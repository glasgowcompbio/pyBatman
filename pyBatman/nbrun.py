# Copyright (c) 2015 Antonino Ingargiola
# License: MIT
# From https://github.com/tritemio/nbrun, slighly modified by Joe

import os
import time

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import HTMLExporter
import codecs

def dict_to_code(mapping):
    lines = ("{} = {}".format(key, repr(value))
             for key, value in mapping.items())
    return '\n'.join(lines)

def run_notebook(nb_name_input, nb_name_output, nb_kwargs=None,
                 insert_pos=1, timeout=3600, execute_kwargs=None):

    timestamp_cell = "**Executed:** %s\n\n**Duration:** %d seconds."
    if nb_kwargs is not None:
        header = '# Cell inserted during automated execution.'
        code = dict_to_code(nb_kwargs)
        code_cell = '\n'.join((header, code))

    if execute_kwargs is None:
        execute_kwargs = {}
    ep = ExecutePreprocessor(timeout=timeout, **execute_kwargs)
    nb = nbformat.read(nb_name_input, as_version=4)
    if len(nb_kwargs) > 0:
        nb['cells'].insert(1, nbformat.v4.new_code_cell(code_cell))

    start_time = time.time()
    try:
        # Execute the notebook
        ep.preprocess(nb, {'metadata': {'path': './'}})
    except:
        # Execution failed, print a message then raise.
        msg = 'Error executing the notebook "%s".\n\n' % nb_name_input
        msg += 'See notebook "%s" for the traceback.' % nb_name_output
        print(msg)
        raise
    else:
        # On successful execution, add timestamping cell
        duration = time.time() - start_time
        timestamp_cell = timestamp_cell % (time.ctime(start_time), duration)
        nb['cells'].insert(0, nbformat.v4.new_markdown_cell(timestamp_cell))
    finally:
        # Save the notebook to HTML output, even when it raises an error
        nbformat.write(nb, nb_name_output)
        exporter = HTMLExporter()
        output, resources = exporter.from_notebook_node(
            nbformat.convert(nb, nbformat.current_nbformat)
        )
        codecs.open(nb_name_output, 'w', encoding='utf-8').write(output)