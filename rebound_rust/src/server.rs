//! server.rs — web server for platform-independent visualization (from
//! server.c; (c) 2023 Hanno Rein, Dave O'Hallaron / Carnegie Mellon).
//!
//! The C serves HTTP from a background thread that dereferences the
//! simulation directly, synchronized with the integrate loop through a
//! shared mutex. Safe Rust cannot alias `&mut reb_simulation` across
//! threads, so the same handshake is expressed with owned shared state:
//! the server thread posts a `need_copy` request and the integrate loop
//! (which holds the simulation) publishes a serialized snapshot and
//! drains the key queue at the same points where the C locks/unlocks
//! its mutex. The HTTP protocol — endpoints, headers, error pages,
//! keyboard commands, screenshot upload — matches server.c.
//!
//! Endpoints: `/` `/index.html` `/rebound.html` (serves rebound.html,
//! auto-downloading it via curl like the C), `/simulation` (binarydata
//! blob of the current state), `/keyboard/<key>`, `/favicon.ico`,
//! `/screenshot` (POST, dataURL base64).
//!
//! Part of rebound_rs, a GPL-3.0-or-later translation of REBOUND 5.1.1.

use crate::binarydata::reb_binarydata_simulation_to_stream;
use crate::tools::{reb_simulation_error, reb_simulation_warning};
use crate::types::*;
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::{Arc, Mutex};

/// server.c `reb_server_header`.
pub const reb_server_header: &str = "HTTP/1.1 200 OK\n\
Server: REBOUND Webserver\n\
Cache-Control: no-cache, no-store, must-revalidate\n\
Pragma: no-cache\n\
Expires: 0\n\
Content-type: text/html\n\
\r\n";

/// server.c `reb_server_header_png`.
pub const reb_server_header_png: &str = "HTTP/1.1 200 OK\n\
Server: REBOUND Webserver\n\
Content-type: image/png\n\
\r\n";

/// rebound.c `reb_favicon_png` (581 bytes).
pub const reb_favicon_png: [u8; 581] = [
    0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d, 0x49, 0x48, 0x44,
    0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10, 0x08, 0x06, 0x00, 0x00, 0x00, 0x1f,
    0xf3, 0xff, 0x61, 0x00, 0x00, 0x00, 0x01, 0x73, 0x52, 0x47, 0x42, 0x00, 0xae, 0xce, 0x1c,
    0xe9, 0x00, 0x00, 0x00, 0x44, 0x65, 0x58, 0x49, 0x66, 0x4d, 0x4d, 0x00, 0x2a, 0x00, 0x00,
    0x00, 0x08, 0x00, 0x01, 0x87, 0x69, 0x00, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
    0x1a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xa0, 0x01, 0x00, 0x03, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x01, 0x00, 0x00, 0xa0, 0x02, 0x00, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
    0x10, 0xa0, 0x03, 0x00, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00,
    0x00, 0x00, 0x34, 0x55, 0x71, 0xf2, 0x00, 0x00, 0x01, 0xaf, 0x49, 0x44, 0x41, 0x54, 0x38,
    0x11, 0x6d, 0xd3, 0x3f, 0x48, 0x55, 0x51, 0x1c, 0xc0, 0xf1, 0xf7, 0xb4, 0x87, 0x29, 0x1a,
    0xa4, 0x90, 0x6d, 0x89, 0xa2, 0x0e, 0x49, 0x69, 0xbd, 0x10, 0x1a, 0x2a, 0x74, 0x14, 0x57,
    0x77, 0x09, 0x5d, 0x0c, 0x23, 0x52, 0x70, 0x11, 0x6c, 0x70, 0x73, 0x30, 0x78, 0x98, 0x36,
    0x94, 0x20, 0x51, 0xad, 0xd2, 0x5f, 0x88, 0x1c, 0x55, 0x5c, 0x24, 0x0c, 0x27, 0x51, 0xb4,
    0xa9, 0x10, 0x41, 0x6a, 0xe8, 0x9f, 0xf8, 0xfd, 0x5e, 0xee, 0x79, 0x1e, 0xd4, 0x1f, 0x7c,
    0xee, 0xf9, 0x77, 0xcf, 0xb9, 0xf7, 0xfc, 0xee, 0xb9, 0xd9, 0xcc, 0xc9, 0xb8, 0x42, 0x57,
    0x1f, 0x6e, 0xe1, 0x07, 0xfe, 0xe2, 0x00, 0x1b, 0x78, 0x82, 0x75, 0x9c, 0x1a, 0x65, 0xf4,
    0x3e, 0x46, 0x01, 0x8d, 0xe8, 0xc5, 0x03, 0x18, 0x25, 0x18, 0xc6, 0x3e, 0x26, 0xe1, 0xbd,
    0x49, 0x38, 0x60, 0x94, 0xe3, 0x43, 0x5a, 0xce, 0x51, 0x6e, 0x63, 0x16, 0xbe, 0xcd, 0x23,
    0xbc, 0x83, 0x71, 0x09, 0xf3, 0x78, 0x81, 0xb3, 0xc8, 0x64, 0xbd, 0x10, 0x53, 0x58, 0xc1,
    0x6f, 0x34, 0xa1, 0x15, 0x7b, 0xd8, 0x41, 0x1b, 0x7c, 0xd0, 0x2f, 0x78, 0xcf, 0x22, 0x9c,
    0xdc, 0x8d, 0xfb, 0x67, 0xb8, 0x5c, 0x83, 0x13, 0x9f, 0x23, 0x8e, 0x09, 0x1a, 0x79, 0x98,
    0x87, 0x87, 0xf8, 0x83, 0x1b, 0xe8, 0x80, 0xf9, 0x69, 0xc1, 0x33, 0xdf, 0x60, 0x06, 0xde,
    0xec, 0x42, 0x9d, 0xf8, 0x87, 0x2d, 0x54, 0x63, 0x04, 0x97, 0xe1, 0xb8, 0x39, 0x58, 0x83,
    0xe1, 0xbc, 0x21, 0x34, 0xd8, 0xf8, 0x8a, 0xcf, 0x30, 0xf3, 0x35, 0xb8, 0x03, 0x33, 0x5f,
    0x8b, 0x10, 0xe7, 0xa9, 0xbc, 0x85, 0xdb, 0x8b, 0xe3, 0xbd, 0x8d, 0x65, 0x54, 0x46, 0xbd,
    0xd3, 0xd4, 0x07, 0x31, 0x16, 0xf5, 0x59, 0xbd, 0x00, 0x13, 0x5d, 0x61, 0x23, 0x8d, 0x79,
    0x93, 0xe3, 0x1e, 0x7f, 0xa6, 0x1d, 0x17, 0x29, 0xcf, 0xa1, 0x80, 0xf6, 0xb4, 0x2f, 0x14,
    0xdf, 0xa9, 0xb8, 0x95, 0xd1, 0xd0, 0x41, 0x99, 0x75, 0x81, 0xd2, 0xa8, 0xa3, 0x87, 0xfa,
    0x4b, 0x78, 0x70, 0x5c, 0xb4, 0x0a, 0x71, 0x7c, 0xa2, 0xe1, 0xbe, 0x7d, 0x90, 0x91, 0x73,
    0x81, 0x6f, 0xa8, 0x83, 0x61, 0x76, 0x17, 0x92, 0x5a, 0x26, 0xf3, 0x85, 0xf2, 0x6a, 0x5a,
    0x8f, 0x0b, 0x93, 0x7e, 0x17, 0xf5, 0xd8, 0x74, 0xc0, 0x6f, 0xee, 0xe9, 0x32, 0x4c, 0x54,
    0x08, 0xbf, 0xf3, 0xbd, 0xd0, 0x88, 0x4a, 0xdf, 0xd8, 0x5c, 0x38, 0xa7, 0xd5, 0x37, 0x58,
    0x45, 0x0e, 0xb7, 0x11, 0x6f, 0xc7, 0x43, 0x93, 0xc7, 0xf1, 0xf8, 0x4f, 0xc7, 0x3e, 0xc2,
    0xdc, 0x64, 0xdc, 0xb3, 0xfd, 0x0a, 0x1f, 0x93, 0xd6, 0xd1, 0xc5, 0x27, 0xc5, 0x8b, 0x3a,
    0xd2, 0x81, 0xd7, 0x48, 0x8e, 0x72, 0x18, 0x74, 0xd5, 0x37, 0xb8, 0x8e, 0x2e, 0xf8, 0xe7,
    0xed, 0xc2, 0x2f, 0x72, 0x13, 0x4b, 0x68, 0xc6, 0x38, 0x1a, 0x30, 0x00, 0x4f, 0x6f, 0xf1,
    0x5f, 0xb0, 0x1e, 0xc2, 0x1f, 0xa8, 0x1f, 0x75, 0xf0, 0xc4, 0xb9, 0x35, 0xcf, 0x8a, 0x07,
    0xee, 0x29, 0xd6, 0x50, 0x8c, 0x43, 0x2d, 0xa4, 0x52, 0x79, 0x8a, 0xe5, 0x13, 0x77, 0x00,
    0x00, 0x00, 0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82,
];
/// rebound.c `reb_favicon_len`.
pub const reb_favicon_len: usize = 581;

/// Shared state between the server thread and the integrate loop.
#[derive(Default)]
pub struct reb_server_shared {
    /// Latest serialized simulation (binarydata format).
    pub snapshot: Vec<u8>,
    /// Incremented whenever `snapshot` is refreshed.
    pub snapshot_seq: u64,
    /// Set by the server thread when a fresh snapshot is wanted
    /// (C: `data->need_copy`).
    pub need_copy: bool,
    /// Keys received via `/keyboard/<n>`, applied by the sim thread.
    pub key_queue: Vec<i32>,
    /// Base64-decoded screenshot received via POST /screenshot.
    pub screenshot: Vec<u8>,
    pub N_screenshot: usize,
    /// C: `data->ready` (1 once listening, -1 on bind error).
    pub ready: i32,
    pub shutdown: bool,
}

/// rebound.h `struct reb_server_data` (thread-handle form).
pub struct reb_server_data {
    pub port: u16,
    pub shared: Arc<Mutex<reb_server_shared>>,
    pub thread: Option<std::thread::JoinHandle<()>>,
}

/// server.c static `base64_decode` (from MIT/FreeBSD wpa utils).
/// Returns None on failure.
pub fn base64_decode(src: &[u8]) -> Option<Vec<u8>> {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut dtable = [0x80u8; 256];
    for (i, &c) in TABLE.iter().enumerate() {
        dtable[c as usize] = i as u8;
    }
    dtable[b'=' as usize] = 0;

    let count = src.iter().filter(|&&c| dtable[c as usize] != 0x80).count();
    if count == 0 || count % 4 != 0 {
        return None;
    }
    let mut out: Vec<u8> = Vec::with_capacity(count / 4 * 3);
    let mut block = [0u8; 4];
    let mut bcount = 0;
    let mut pad = 0;
    for &c in src {
        let tmp = dtable[c as usize];
        if tmp == 0x80 {
            continue;
        }
        if c == b'=' {
            pad += 1;
        }
        block[bcount] = tmp;
        bcount += 1;
        if bcount == 4 {
            out.push((block[0] << 2) | (block[1] >> 4));
            out.push((block[1] << 4) | (block[2] >> 2));
            out.push((block[2] << 6) | block[3]);
            bcount = 0;
            if pad > 0 {
                if pad == 1 {
                    out.pop();
                } else if pad == 2 {
                    out.pop();
                    out.pop();
                } else {
                    /* Invalid padding */
                    return None;
                }
                break;
            }
        }
    }
    Some(out)
}

/// server.c static `reb_server_cerror`.
fn reb_server_cerror(stream: &mut TcpStream, cause: &str) {
    let buf = format!(
        "HTTP/1.1 501 Not Implemented\n\
Content-type: text/html\n\
\n\
<html><title>REBOUND Webserver Error</title><body>\n\
<h1>Error</h1>\n\
<p>{}</p>\n\
<hr><em>REBOUND Webserver</em>\n\
</body></html>\n",
        cause
    );
    println!("\nREBOUND Webserver error: {}", cause);
    let _ = stream.write_all(buf.as_bytes());
}

/// Wait (briefly) for the integrate loop to publish a fresh snapshot;
/// returns whatever snapshot is available afterwards. This is the safe
/// Rust rendering of the C's need_copy/mutex handshake.
fn request_snapshot(shared: &Arc<Mutex<reb_server_shared>>) -> Vec<u8> {
    let seq0;
    {
        let mut s = shared.lock().unwrap();
        seq0 = s.snapshot_seq;
        s.need_copy = true;
    }
    for _ in 0..200 {
        {
            let s = shared.lock().unwrap();
            if s.snapshot_seq != seq0 || s.shutdown {
                return s.snapshot.clone();
            }
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    shared.lock().unwrap().snapshot.clone()
}

/// server.c `reb_server_start` — the server thread main loop.
fn reb_server_thread(listener: TcpListener, shared: Arc<Mutex<reb_server_shared>>, port: u16) {
    println!(
        "REBOUND Webserver listening on http://localhost:{} (not secure) ...",
        port
    );
    for conn in listener.incoming() {
        if shared.lock().unwrap().shutdown {
            break;
        }
        let mut stream = match conn {
            Ok(s) => s,
            Err(_) => break, // Accept fails when the socket is closed.
        };
        // Receive request (headers + optional body).
        let mut recbuf: Vec<u8> = Vec::new();
        let mut chunk = [0u8; 1024];
        loop {
            match stream.read(&mut chunk) {
                Ok(0) => break,
                Ok(n) => {
                    recbuf.extend_from_slice(&chunk[..n]);
                    if n < chunk.len() {
                        break;
                    }
                }
                Err(_) => break,
            }
        }
        let text = String::from_utf8_lossy(&recbuf).to_string();
        let mut lines = text.split('\n');
        let request_line = lines.next().unwrap_or("");
        let mut parts = request_line.split_whitespace();
        let method = parts.next().unwrap_or("");
        let uri = parts.next().unwrap_or("");

        // only support the GET and POST methods
        if !method.eq_ignore_ascii_case("GET") && !method.eq_ignore_ascii_case("POST") {
            reb_server_cerror(&mut stream, "Only GET+POST are implemented.");
            continue;
        }

        // read (and otherwise ignore) the HTTP headers
        let mut content_length: usize = 0;
        let mut body_start = 0usize;
        if let Some(p) = text.find("\r\n\r\n") {
            body_start = p + 4;
        }
        for line in text[..body_start.min(text.len())].lines() {
            if let Some(cl) = line.strip_prefix("Content-Length: ") {
                content_length = cl.trim().parse().unwrap_or(0);
            }
        }

        if uri.eq_ignore_ascii_case("/simulation") {
            let snapshot = request_snapshot(&shared);
            let _ = stream.write_all(reb_server_header.as_bytes());
            let _ = stream.write_all(&snapshot);
        } else if uri.len() > 10 && uri[..10].eq_ignore_ascii_case("/keyboard/") {
            let key: i32 = uri[10..].parse().unwrap_or(0);
            shared.lock().unwrap().key_queue.push(key);
            let _ = stream.write_all(reb_server_header.as_bytes());
            let _ = stream.write_all(b"ok.\n");
        } else if uri.eq_ignore_ascii_case("/")
            || uri.eq_ignore_ascii_case("/index.html")
            || uri.eq_ignore_ascii_case("/rebound.html")
        {
            match std::fs::read("rebound.html") {
                Ok(html) => {
                    let _ = stream.write_all(reb_server_header.as_bytes());
                    let _ = stream.write_all(&html);
                }
                Err(_) => {
                    reb_server_cerror(
                        &mut stream,
                        "rebound.html not found in current directory. Try `make rebound.html`.",
                    );
                }
            }
        } else if uri.eq_ignore_ascii_case("/favicon.ico") {
            let _ = stream.write_all(reb_server_header_png.as_bytes());
            let _ = stream.write_all(&reb_favicon_png);
        } else if uri.eq_ignore_ascii_case("/screenshot") {
            if content_length == 0 {
                println!("Received screenshot with size zero.");
            } else {
                let body = &recbuf[body_start.min(recbuf.len())..];
                let data_url = String::from_utf8_lossy(body);
                match data_url.find(',') {
                    None => {
                        println!("Unable to decode received screenshot. Data not in dataURL format.")
                    }
                    Some(cpos) => {
                        let b64 = &data_url[cpos + 1..];
                        match base64_decode(b64.trim_end_matches('\0').as_bytes()) {
                            Some(png) => {
                                let mut s = shared.lock().unwrap();
                                s.N_screenshot = png.len();
                                s.screenshot = png;
                                s.key_queue.push(-1); // internal: pause after screenshot
                            }
                            None => println!("An error occurred while decoding the screenshot."),
                        }
                    }
                }
            }
            let _ = stream.write_all(reb_server_header.as_bytes());
            let _ = stream.write_all(b"ok.\n");
        } else {
            reb_server_cerror(&mut stream, "Unsupported URI.");
            println!("URI: {}", uri);
        }
        let _ = stream.flush();
    }
    println!("Server shutting down...");
}

/// Called from the integrate loop (and from the pause wait loop):
/// publishes a snapshot when the server requested one and applies
/// queued keyboard commands (the C server thread applies these directly
/// under the shared mutex).
pub fn reb_server_update(r: &mut reb_simulation) {
    let shared = match &r.server_data {
        Some(sd) => sd.shared.clone(),
        None => return,
    };
    let (need_copy, keys) = {
        let mut s = shared.lock().unwrap();
        (s.need_copy, std::mem::take(&mut s.key_queue))
    };
    for key in keys {
        let skip_default_keys = 0; // key_callback (SERVER/OPENGL builds) not carried
        if skip_default_keys == 0 {
            match key {
                81 => {
                    // 'Q'
                    r.status = REB_STATUS_USER;
                }
                32 => {
                    // ' '
                    if r.status == REB_STATUS_PAUSED {
                        println!("Resume.");
                        r.status = REB_STATUS_RUNNING;
                    } else if r.status == REB_STATUS_RUNNING {
                        println!("Pause.");
                        r.status = REB_STATUS_PAUSED;
                    }
                }
                264 => {
                    // arrow down
                    if r.status == REB_STATUS_PAUSED {
                        r.status = REB_STATUS_SINGLE_STEP;
                        println!("Step.");
                    }
                }
                267 => {
                    // page down
                    if r.status == REB_STATUS_PAUSED {
                        r.status = REB_STATUS_SINGLE_STEP - 50;
                        println!("50 steps.");
                    }
                }
                -1 => {
                    // internal: screenshot received
                    r.status = REB_STATUS_PAUSED;
                }
                _ => {}
            }
        }
    }
    if need_copy {
        let snapshot = reb_binarydata_simulation_to_stream(r);
        let mut s = shared.lock().unwrap();
        s.snapshot = snapshot;
        s.snapshot_seq += 1;
        s.need_copy = false;
    }
}

/// server.c `reb_simulation_start_server`. Returns 0 on success, -1 on
/// error (same messages as the C).
pub fn reb_simulation_start_server(r: &mut reb_simulation, port: i32) -> i32 {
    if port == 0 {
        reb_simulation_error(r, "Cannot start server. Invalid port.");
        return -1;
    }
    if r.server_data.is_some() {
        reb_simulation_error(r, "Server already started.");
        return -1;
    }
    if !std::path::Path::new("rebound.html").exists() {
        reb_simulation_warning(
            r,
            "File rebound.html not found in current directory. Attempting to download it from github.",
        );
        let _ = std::process::Command::new("curl")
            .args([
                "-L",
                "-s",
                "--output",
                "rebound.html",
                "https://github.com/hannorein/rebound/releases/latest/download/rebound.html",
            ])
            .status();
        if !std::path::Path::new("rebound.html").exists() {
            reb_simulation_warning(r, "Automatic download failed. Manually download the file from github and place it in the current directory to enable browser based visualization.");
        } else {
            println!("Success: rebound.html downloaded.");
        }
    }
    let listener = match TcpListener::bind(("0.0.0.0", port as u16)) {
        Ok(l) => l,
        Err(_) => {
            let msg = format!("Error binding to port {}. Port might be in use.\n", port);
            reb_simulation_error(r, &msg);
            return -1;
        }
    };
    let shared = Arc::new(Mutex::new(reb_server_shared::default()));
    shared.lock().unwrap().ready = 1;
    let shared_thread = shared.clone();
    let thread = std::thread::spawn(move || {
        reb_server_thread(listener, shared_thread, port as u16);
    });
    r.server_data = Some(reb_server_data {
        port: port as u16,
        shared,
        thread: Some(thread),
    });
    0
}

/// server.c `reb_simulation_stop_server`.
pub fn reb_simulation_stop_server(r: &mut reb_simulation) {
    if let Some(mut sd) = r.server_data.take() {
        sd.shared.lock().unwrap().shutdown = true;
        // Unblock accept() by connecting once (the C closes the socket).
        let _ = TcpStream::connect(("127.0.0.1", sd.port));
        if let Some(t) = sd.thread.take() {
            let _ = t.join();
        }
    }
}
