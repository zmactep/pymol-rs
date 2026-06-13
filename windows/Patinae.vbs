Set WshShell = CreateObject("WScript.Shell")
Set FSO = CreateObject("Scripting.FileSystemObject")

BundleDir = FSO.GetParentFolderName(WScript.ScriptFullName) & "\"

WshShell.Environment("Process")("PATINAE_PLUGIN_DIR") = BundleDir & "plugins"
WshShell.Environment("Process")("VIRTUAL_ENV") = BundleDir & "python-venv"
WshShell.Environment("Process")("PATH") = BundleDir & "python;" & WshShell.Environment("Process")("PATH")

Dim Args
Args = ""
For Each a In WScript.Arguments
    Args = Args & " """ & a & """"
Next

WshShell.Run """" & BundleDir & "patinae.exe""" & Args, 0, False
