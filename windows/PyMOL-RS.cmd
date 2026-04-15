@echo off
set "BUNDLE=%~dp0"
set "PYMOL_RS_PLUGIN_DIR=%BUNDLE%plugins"
set "VIRTUAL_ENV=%BUNDLE%python-venv"
set "PATH=%BUNDLE%python;%PATH%"
start "" "%BUNDLE%pymol-rs.exe" %*
