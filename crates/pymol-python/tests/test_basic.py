"""
Basic tests for pymol-python bindings.

These tests verify the core functionality of the Python bindings.
"""

import pytest


def test_import():
    """Test that the module can be imported."""
    import pymol
    from pymol import cmd
    assert cmd is not None


def test_version():
    """Test version information."""
    from pymol._pymol import __version__
    assert isinstance(__version__, str)
    assert len(__version__) > 0


def test_cmd_exists():
    """Test that cmd object has expected methods."""
    from pymol import cmd
    
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
    from pymol import Color
    
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
    from pymol import Element
    
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
    from pymol import SelectionResult
    
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
    from pymol import ObjectMolecule
    
    mol = ObjectMolecule("test")
    assert mol.name == "test"
    assert mol.atom_count == 0
    assert mol.bond_count == 0


def test_io_module():
    """Test io module functions exist."""
    from pymol import io
    
    assert hasattr(io, 'read_file')
    assert hasattr(io, 'write_file')
    assert hasattr(io, 'read_all')
    assert hasattr(io, 'fetch')
    assert hasattr(io, 'detect_format')


def test_selecting_module():
    """Test selecting module functions exist."""
    from pymol import selecting
    
    assert hasattr(selecting, 'select')
    assert hasattr(selecting, 'select_atoms')
    assert hasattr(selecting, 'count_atoms')


def test_color_module():
    """Test color module functionality."""
    from pymol import color
    
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
    from pymol import settings
    
    # List settings
    all_settings = settings.list_settings()
    assert isinstance(all_settings, list)
    assert len(all_settings) > 0
    
    # Search settings
    cartoon_settings = settings.search_settings("cartoon")
    assert isinstance(cartoon_settings, list)


def test_exceptions():
    """Test exception types exist."""
    from pymol import PymolError, SelectionError
    
    assert PymolError is not None
    assert SelectionError is not None
    
    # SelectionError should be a subclass of PymolError
    assert issubclass(SelectionError, PymolError)


def test_png_command_headless():
    """Test that png command exists and raises appropriate error in headless mode."""
    from pymol import cmd
    
    # In headless mode, png should raise an error explaining the limitation
    with pytest.raises(RuntimeError, match="GUI mode"):
        cmd.png("test_output.png")


def test_png_with_dimensions():
    """Test png command signature accepts width and height."""
    from pymol import cmd
    
    # Should accept width and height parameters
    # In headless mode this will raise, but we're testing the API exists
    with pytest.raises(RuntimeError):
        cmd.png("test.png", width=800, height=600)
    
    with pytest.raises(RuntimeError):
        cmd.png("test.png", 1920, 1080)


def test_show_viewer_exists():
    """Test that show_viewer method exists and has correct signature."""
    from pymol import cmd
    
    # Check method exists
    assert hasattr(cmd, 'show_viewer')
    
    # Check it's callable
    assert callable(cmd.show_viewer)


def test_viewer_helper_methods_exist():
    """Test that viewer helper methods exist."""
    from pymol import cmd
    
    # Check methods exist
    assert hasattr(cmd, 'viewer_is_open')
    assert hasattr(cmd, 'close_viewer')
    assert hasattr(cmd, 'wait_viewer')
    
    # Check they're callable
    assert callable(cmd.viewer_is_open)
    assert callable(cmd.close_viewer)
    assert callable(cmd.wait_viewer)


def test_viewer_not_running_initially():
    """Test that viewer_is_open returns False initially."""
    from pymol import cmd
    
    # Viewer should not be running initially
    assert cmd.viewer_is_open() == False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
