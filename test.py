from pymol_rs import cmd

def load_and_prepare(structure):
    cmd.fetch(structure, sync=True)
    cmd.show_as("cartoon")
    cmd.color("by_chain")

cmd.extend("load_and_prepare", load_and_prepare)
