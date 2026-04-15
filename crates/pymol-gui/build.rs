fn main() {
    // Embed icon and version info in the Windows executable
    #[cfg(target_os = "windows")]
    {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("../../images/pymol-rs.ico");
        res.set("ProductName", "PyMOL-RS");
        res.set("FileDescription", "PyMOL-RS Molecular Visualization");
        res.compile().expect("Failed to compile Windows resources");
    }
}
