use sha2::{Digest, Sha256};
fn main() {
    let rustc = std::env::var_os("RUSTC").expect("Cargo supplies RUSTC");
    let output = std::process::Command::new(rustc)
        .arg("--version")
        .output()
        .expect("rustc version");
    assert!(output.status.success());
    println!(
        "cargo:rustc-env=CRN_RUSTC_VERSION={}",
        String::from_utf8(output.stdout)
            .expect("UTF-8 rustc version")
            .trim()
    );
    // Normalized crates have their own lockfile; a checkout uses the workspace.
    let root = std::path::PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let local = root.join("Cargo.lock");
    let lock = if local.is_file() {
        local
    } else {
        root.join("../../Cargo.lock")
    };
    println!("cargo:rerun-if-changed={}", lock.display());
    let bytes = std::fs::read(lock).expect("committed or packaged Cargo.lock");
    println!(
        "cargo:rustc-env=CRN_CARGO_LOCK_SHA256={:x}",
        Sha256::digest(bytes)
    );
    println!("cargo:rerun-if-changed=build.rs");
}
