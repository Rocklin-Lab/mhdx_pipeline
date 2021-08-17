"""Example Google style docstrings.

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html

"""
import os
import shutil
import argparse

def main(input_paths, output_paths):
	"""Moves idotp_check passing tensors to a new directory to simplify downstream input calling.

    Args:
        input_paths (list of strings): Sorted list of paths/to/inputs in the same order as the output list.
        output_paths (list of strings): Sorted list of paths/to/outputs in the same order as the input list.

    Returns:
        None

    """
	for fin, fout in zip(input_paths, output_paths):
		os.makedirs(os.path.dirname(fout), exist_ok=True)
		shutil.copy(fin, fout)

if __name__ == "__main__":

	if "snakemake" in globals():
		main(snakemake.input, snakemake.output)

	else:
		parser = argparse.ArgumentParser()
		parser.add_argument("library_info_path", help="path/to/checked_library_info.json")
		parser.add_argument("input_dir_path", help="path/to/input_dir, inputs are globbed from this path, I don't work right now", default="resources/")
		parser.add_argument("output_dir_path", help="path/to/output_dir, outputs are globbed from this path, I don't work right now")
		args = parser.parse_args()

		# Make explicit inputs and outputs from checked_library_info.json.
