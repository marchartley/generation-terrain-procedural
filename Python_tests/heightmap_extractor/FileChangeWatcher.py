import os
import sys
import time
import multiprocessing


class FileChangeWatcher(object):
    running = True
    refresh_delay_secs = 1
    _thread = None

    # Constructor
    def __init__(self, watch_file, call_func_on_change=None, *args, **kwargs):
        self._cached_stamp = 0
        self.filename = watch_file
        self.call_func_on_change = call_func_on_change
        self.args = args
        self.kwargs = kwargs

    # Look for changes
    def look(self):
        stamp = os.stat(self.filename).st_mtime
        if stamp != self._cached_stamp:
            self._cached_stamp = stamp
            # File has changed, so do something...
            if self.call_func_on_change is not None:
                self.call_func_on_change(*self.args, **self.kwargs)

    # Keep watching in a loop
    def watch(self):
        while self.running:
            try:
                # Look for changes
                time.sleep(self.refresh_delay_secs)
                self.look()
            except KeyboardInterrupt:
                break
            except FileNotFoundError:
                # Action on file not found
                pass
            # except:
            #     print('Unhandled error: %s' % sys.exc_info()[0])

    def start(self):
        self.running = True
        self._thread = multiprocessing.Process(target = self.watch)
        self._thread.start()

    def stop(self):
        self.running = False
        self._thread.terminate()


if __name__ == "__main__":
    def testFunction(text_to_print):
        print(text_to_print)

    watcher = FileChangeWatcher("test_fourier_image.png", testFunction, text_to_print="File changed!")
    watcher.start()
    time.sleep(60)
    watcher.stop()
