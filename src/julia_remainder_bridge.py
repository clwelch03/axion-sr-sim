"""Experimental single-caller Julia bridge; not enabled by the solver.

One persistent subprocess, two shared float64 buffers, one input copy and one
output copy per call. Returned arrays own their data. No third-party bridge
packages or Julia process launches per RHS. Use as a context manager.
"""

from pathlib import Path
import select
import subprocess
import tempfile

import numpy as np


class JuliaRemainder:
    def __init__(self, julia, capacity=2048 * 256):
        if not isinstance(capacity, int) or capacity <= 0:
            raise ValueError("capacity must be a positive integer")
        self.capacity = capacity
        self.directory = tempfile.TemporaryDirectory(prefix='axion-julia-bridge-')
        self.process = None
        self.input = self.output = None
        try:
            input_path = Path(self.directory.name) / 'input.bin'
            output_path = Path(self.directory.name) / 'output.bin'
            self.input = np.memmap(input_path, mode='w+', dtype='float64', shape=(capacity,))
            self.output = np.memmap(output_path, mode='w+', dtype='float64', shape=(capacity,))
            self.process = subprocess.Popen(
                [julia, '--startup-file=no', '--history-file=no', '--threads=1',
                 str(Path(__file__).with_name('julia_remainder_worker.jl')),
                 str(input_path), str(output_path), str(capacity)],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
            self._expect('ready', 120)
        except BaseException:
            self.close()
            raise

    def _expect(self, expected, timeout):
        if not select.select([self.process.stdout], [], [], timeout)[0]:
            raise TimeoutError('Julia worker did not respond')
        response = self.process.stdout.readline().strip()
        if response != expected:
            raise RuntimeError(f'Julia worker returned {response!r}, expected {expected!r}')

    def __call__(self, field):
        if self.process is None:
            raise RuntimeError('Julia worker is closed')
        field = np.asarray(field)
        if np.iscomplexobj(field):
            raise ValueError('the axion potential requires a real field')
        field = np.asarray(field, dtype=float)
        if field.size > self.capacity:
            raise ValueError('field exceeds shared buffer capacity')
        # Flattening is purely storage order: this is an elementwise function.
        self.input[:field.size] = field.ravel()
        self.process.stdin.write(f'{field.size}\n')
        self.process.stdin.flush()
        self._expect('done', 30)
        return np.array(self.output[:field.size]).reshape(field.shape)

    def close(self):
        if self.process is not None:
            if self.process.poll() is None:
                try:
                    self.process.stdin.write('quit\n')
                    self.process.stdin.flush()
                except BrokenPipeError:
                    pass
                try:
                    self.process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    self.process.wait()
            self.process.stdin.close()
            self.process.stdout.close()
            self.process = None
        for buffer in (self.input, self.output):
            if buffer is not None:
                buffer._mmap.close()
        self.input = self.output = None
        self.directory.cleanup()

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.close()
