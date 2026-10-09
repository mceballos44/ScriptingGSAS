# Keeps GSAS-II's printouts out of the terminal.
#
# Everything printed inside gsas_output_to() goes to a log file. Lines that
# look like problems (contain one of the keywords) are also shown in the
# terminal, prefixed with [GSAS], so warnings and failures are not missed.

### Written by Mauricio Ceballos, Joester Group, Northwestern University

import sys
import contextlib


class _LogStream:
    """
    Stand-in for sys.stdout: writes everything to the log file and passes
    matching lines through to the terminal
    """
    def __init__(self, log, terminal, keywords):
        self.log = log
        self.terminal = terminal
        self.keywords = [k.lower() for k in keywords]
        self._partial = ''

    def write(self, text):
        self.log.write(text)
        # Only check complete lines; keep the unfinished end for next time
        self._partial += text
        *lines, self._partial = self._partial.split('\n')
        for line in lines:
            self._show_if_important(line)
        return len(text)

    def _show_if_important(self, line):
        if any(k in line.lower() for k in self.keywords):
            self.terminal.write(f'    [GSAS] {line.strip()}\n')

    def finish(self):
        if self._partial:
            self._show_if_important(self._partial)
            self._partial = ''
        self.flush()

    def flush(self):
        self.log.flush()
        self.terminal.flush()


@contextlib.contextmanager
def gsas_output_to(log_file, keywords):
    """
    Send everything printed inside the with-block to log_file.

    Inputs:
    log_file: Path of the log (overwritten, folder created if missing)
    keywords: lines containing any of these (case-insensitive) are also
        shown in the terminal

    Usage:
        with gsas_output_to(LOG_DIR / 'AFP3.log', GSAS_LOG_KEYWORDS):
            full_analysis('AFP3')

    Errors (exceptions) still reach the terminal. Output written directly
    by GSAS-II's compiled code (not through Python) is not captured.
    """
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, 'w') as log:
        stream = _LogStream(log, sys.stdout, keywords)
        try:
            with contextlib.redirect_stdout(stream):
                yield
        finally:
            stream.finish()
