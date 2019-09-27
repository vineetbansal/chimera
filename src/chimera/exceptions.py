import sys
import platform
import struct
import subprocess
import traceback


def handle_exception(exc_type, exc_value, exc_traceback):
    """
    Handle any top-level unhandled exception.
    :param exc_type: Exception type object
    :param exc_value: Exception value object (an instance of type exc_type)
    :param exc_traceback: The Traceback object associated with exc_value (also available as exc_value.__traceback__)
    :return: On return, useful diagnostic information has been logged, and the exception re-raised.
    """

    # Are we explicitly/interactively killing a run? Just do it.
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    from chimera.utils import get_full_version

    lines = list()

    lines.append(f'Application version: {get_full_version()}')
    lines.append(f'Platform: {platform.platform()}')
    lines.append(f'Python version: {sys.version}')
    lines.append(f'Python 32/64 bit: {8 * struct.calcsize("P")}')

    lines.append('conda list output:')
    try:
        lines.extend(subprocess.check_output(['conda', 'list'], stderr=subprocess.STDOUT).decode('utf8').split('\n'))
    except:  # nopep8
        pass

    lines.append('pip freeze output:')
    try:
        lines.extend(subprocess.check_output(['pip', 'freeze'], stderr=subprocess.STDOUT).decode('utf8').split('\n'))
    except:  # nopep8
        pass

    # Walk through all traceback objects (oldest call -> most recent call), capturing frame/local variable information.
    lines.append('Exception Details (most recent call last)')
    frame_generator = traceback.walk_tb(exc_traceback)
    stack_summary = traceback.StackSummary.extract(frame_generator, capture_locals=True)
    frame_strings = stack_summary.format()
    for s in frame_strings:
        lines.extend(s.split('\n'))

    try:
        with open('app.err.log', 'w') as f:
            f.write('\n'.join(lines))
    except:  # nopep8
        pass

    try:
        # re-raise the exception we got for the caller.
        raise exc_value
    finally:
        # cleanup - see https://cosmicpercolator.com/2016/01/13/exception-leaks-in-python-2-and-3/
        del exc_value, exc_traceback