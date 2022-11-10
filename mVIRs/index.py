import logging
import os
import pathlib
import shutil
import subprocess


def index_genome(seq_file: str, output_folder: str) -> None:
    logging.info(f'Start building bwa index on {seq_file}')
    # make output folder
    out_folder = pathlib.Path(output_folder)
    out_folder.mkdir(parents=True, exist_ok=True)
    shutil.copy2(seq_file, output_folder)
    index_path = os.path.abspath(os.path.join(output_folder, seq_file.split("/")[-1]))
    command = f'bwa index {seq_file} -p {index_path}'

    try:
        returncode: int = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))

    logging.info(f'Successfully built index on {seq_file}')