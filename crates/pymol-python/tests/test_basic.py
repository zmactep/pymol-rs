"""
Basic tests for pymol-rs Python bindings.

These tests verify the core functionality of the Python bindings.
"""

import pytest


def test_import():
    """Test that the module can be imported."""
    import pymol_rs
    from pymol_rs import cmd
    assert cmd is not None


def test_version():
    """Test version information."""
    from pymol_rs._pymol_rs import __version__
    assert isinstance(__version__, str)
    assert len(__version__) > 0


def test_cmd_exists():
    """Test that cmd object has expected methods."""
    from pymol_rs import cmd
    
    # Check core methods exist
    assert hasattr(cmd, 'load')
    assert hasattr(cmd, 'save')
    assert hasattr(cmd, 'fetch')
    assert hasattr(cmd, 'show')
    assert hasattr(cmd, 'hide')
    assert hasattr(cmd, 'color')
    assert hasattr(cmd, 'select')
    assert hasattr(cmd, 'zoom')
    assert hasattr(cmd, 'center')
    assert hasattr(cmd, 'reset')
    assert hasattr(cmd, 'delete')
    assert hasattr(cmd, 'get_names')
    assert hasattr(cmd, 'png')
    assert hasattr(cmd, 'show_viewer')
    assert hasattr(cmd, 'viewer_is_open')
    assert hasattr(cmd, 'close_viewer')
    assert hasattr(cmd, 'wait_viewer')


def test_color_class():
    """Test Color class functionality."""
    from pymol_rs import Color
    
    # Create color
    c = Color(1.0, 0.0, 0.0)
    assert c.r == 1.0
    assert c.g == 0.0
    assert c.b == 0.0
    
    # From RGB8
    c2 = Color.from_rgb8(255, 128, 0)
    assert c2.r == pytest.approx(1.0, abs=0.01)
    assert c2.g == pytest.approx(0.5, abs=0.01)
    assert c2.b == 0.0
    
    # From hex
    c3 = Color.from_hex("#ff0000")
    assert c3.r == 1.0


def test_element_class():
    """Test Element class functionality."""
    from pymol_rs import Element
    
    # Create from atomic number
    carbon = Element(6)
    assert carbon.symbol == "C"
    assert carbon.name == "Carbon"
    assert carbon.atomic_number == 6
    
    # Create from symbol
    oxygen = Element.from_symbol("O")
    assert oxygen is not None
    assert oxygen.atomic_number == 8


def test_selection_result():
    """Test SelectionResult class."""
    from pymol_rs import SelectionResult
    
    # Create empty selection
    sel = SelectionResult.none(100)
    assert sel.count == 0
    assert sel.is_empty()
    
    # Create full selection
    sel2 = SelectionResult.all(100)
    assert sel2.count == 100
    assert not sel2.is_empty()
    
    # Create from indices
    sel3 = SelectionResult.from_indices([0, 1, 2, 3], 100)
    assert sel3.count == 4
    assert sel3.contains(0)
    assert sel3.contains(3)
    assert not sel3.contains(50)


def test_object_molecule():
    """Test ObjectMolecule class."""
    from pymol_rs import ObjectMolecule
    
    mol = ObjectMolecule("test")
    assert mol.name == "test"
    assert mol.atom_count == 0
    assert mol.bond_count == 0


def test_io_module():
    """Test io module functions exist."""
    from pymol_rs import io
    
    assert hasattr(io, 'read_file')
    assert hasattr(io, 'write_file')
    assert hasattr(io, 'read_all')
    assert hasattr(io, 'fetch')
    assert hasattr(io, 'detect_format')


def test_selecting_module():
    """Test selecting module functions exist."""
    from pymol_rs import selecting
    
    assert hasattr(selecting, 'select')
    assert hasattr(selecting, 'select_atoms')
    assert hasattr(selecting, 'count_atoms')


def test_color_module():
    """Test color module functionality."""
    from pymol_rs import color
    
    assert hasattr(color, 'get_color')
    assert hasattr(color, 'element_color')
    assert hasattr(color, 'chain_color')
    
    # Check color constants
    assert color.RED is not None
    assert color.GREEN is not None
    assert color.BLUE is not None
    assert color.WHITE is not None
    assert color.BLACK is not None


def test_settings_module():
    """Test settings module."""
    from pymol_rs import settings
    
    # List settings
    all_settings = settings.list_settings()
    assert isinstance(all_settings, list)
    assert len(all_settings) > 0
    
    # Search settings
    cartoon_settings = settings.search_settings("cartoon")
    assert isinstance(cartoon_settings, list)


def test_exceptions():
    """Test exception types exist."""
    from pymol_rs import PymolError, SelectionError
    
    assert PymolError is not None
    assert SelectionError is not None
    
    # SelectionError should be a subclass of PymolError
    assert issubclass(SelectionError, PymolError)


def test_png_command_headless():
    """Test that png command works in headless mode."""
    from pymol_rs import cmd
    import tempfile
    import os
    
    # PNG should now work in headless mode (hidden window with render context)
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "test_output.png")
        cmd.png(output_path)  # Should succeed
        assert os.path.exists(output_path), "PNG file was not created"
        # Verify it's a valid PNG (starts with PNG magic bytes)
        with open(output_path, 'rb') as f:
            magic = f.read(8)
            assert magic[:4] == b'\x89PNG', "File is not a valid PNG"


def test_png_with_dimensions():
    """Test png command with custom width and height."""
    from pymol_rs import cmd
    import tempfile
    import os
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test with width and height keyword arguments
        output_path = os.path.join(tmpdir, "test_800x600.png")
        cmd.png(output_path, width=800, height=600)
        assert os.path.exists(output_path), "PNG file was not created"
        
        # Test with positional arguments
        output_path2 = os.path.join(tmpdir, "test_1920x1080.png")
        cmd.png(output_path2, 1920, 1080)
        assert os.path.exists(output_path2), "PNG file was not created"


def test_show_gui_exists():
    """Test that show_gui method exists and has correct signature."""
    from pymol_rs import cmd
    
    # Check method exists
    assert hasattr(cmd, 'show_gui')
    
    # Check it's callable
    assert callable(cmd.show_gui)


def test_gui_helper_methods_exist():
    """Test that GUI helper methods exist."""
    from pymol_rs import cmd
    
    # Check methods exist
    assert hasattr(cmd, 'gui_is_open')
    assert hasattr(cmd, 'close_gui')
    assert hasattr(cmd, 'wait_gui')
    
    # Check they're callable
    assert callable(cmd.gui_is_open)
    assert callable(cmd.close_gui)
    assert callable(cmd.wait_gui)


def test_gui_not_running_initially():
    """Test that gui_is_open returns False initially."""
    from pymol_rs import cmd
    
    # GUI should not be running initially
    assert cmd.gui_is_open() == False


def test_extend_method_exists():
    """Test that extend method exists for registering custom commands."""
    from pymol_rs import cmd
    
    assert hasattr(cmd, 'extend')
    assert callable(cmd.extend)


def test_unextend_method_exists():
    """Test that unextend method exists for unregistering custom commands."""
    from pymol_rs import cmd
    
    assert hasattr(cmd, 'unextend')
    assert callable(cmd.unextend)


def test_run_method_exists():
    """Test that run method exists for script execution."""
    from pymol_rs import cmd
    
    assert hasattr(cmd, 'run')
    assert callable(cmd.run)


def test_do_method_exists():
    """Test that do method exists for command execution."""
    from pymol_rs import cmd
    
    # The method is named 'do' in Python (with underscore removed)
    assert hasattr(cmd, 'do')
    assert callable(getattr(cmd, 'do'))


def test_extend_custom_command():
    """Test registering a custom command with extend."""
    from pymol_rs import cmd
    
    # Define a simple command
    def my_test_command(arg1="default"):
        return f"Got: {arg1}"
    
    # Register it - should not raise
    cmd.extend("test_cmd", my_test_command)
    
    # Unregister it - should not raise
    cmd.unextend("test_cmd")


def test_do_command_basic():
    """Test executing basic commands via do."""
    from pymol_rs import cmd
    
    # These commands should work in headless mode
    cmd.do("reset")  # Should not raise
    
    # Unknown commands should raise
    with pytest.raises(ValueError, match="Unknown command"):
        cmd.do("nonexistent_command arg1 arg2")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
