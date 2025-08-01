from tracemalloc import Traceback

bl_info = {
    "name": "Auto Reload External Python Files",
    "author": "Marc Hartley",
    "version": (1, 1),
    "blender": (2, 80, 0),
    "description": "Automatically reloads external Python scripts when changed on disk.",
    "category": "Development",
}

import bpy
import os
import subprocess
import importlib
import traceback
from bpy.app.handlers import persistent

_observers = []
_watched_files = set()


def logger(content):
    print(f"[Auto-Reload] {content}")


def ensure_watchdog_installed():
    try:
        import watchdog
        return True
    except ImportError:
        logger("watchdog not found. Attempting installation...")
        try:
            subprocess.check_call([bpy.app.binary_path_python, "-m", "pip", "install", "watchdog"])
            importlib.invalidate_caches()
            import watchdog
            logger("Successfully installed watchdog.")
            return True
        except Exception as e:
            logger("Failed to install watchdog:", e)
            show_warning("Could not install 'watchdog'. See console for details.")
            return False


def show_warning(text):
    def draw(self, context):
        self.layout.label(text=text)

    try:
        bpy.context.window_manager.popup_menu(draw, title="Auto-Reload Add-on", icon='ERROR')
    except Exception as e:
        logger("Didn't manage to show popup")


def find_external_python_files():
    external_files = set()
    logger(f"Checking all files in {bpy.data.texts}...")
    for text in bpy.data.texts:
        logger(f"Checking {text} ?")
        if text.filepath and os.path.isfile(bpy.path.abspath(text.filepath)):
            full_path = bpy.path.abspath(text.filepath)
            if full_path.endswith(".py"):
                logger(f"Added {text} ({full_path})!")
                external_files.add(full_path)
    logger(f"Found {len(external_files)} files to watch")
    return sorted(external_files)


def run_script_later(text_block):
    def _exec_later():
        try:
            win = bpy.context.window_manager.windows[0]
            area = next(area for area in win.screen.areas if area.type == 'TEXT_EDITOR')

            override = {
                'window': win,
                'screen': win.screen,
                'area': area,
                'region': next(r for r in area.regions if r.type == 'WINDOW'),
                'edit_text': text_block,
            }

            bpy.ops.text.run_script(override)
            logger(f"Executed script: {text_block.name}")
        except Exception as e:
            logger(f"Script execution failed: {e}")
        return None  # Don't repeat

    bpy.app.timers.register(_exec_later)


def start_file_watchers():
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler

    class ReloadOnChangeHandler(FileSystemEventHandler):
        def __init__(self, file_path):
            self.file_path = os.path.abspath(file_path)

        def on_modified(self, event):
            if os.path.abspath(event.src_path) == self.file_path:
                logger(f"Detected change in: {self.file_path}")
                try:
                    for text in bpy.data.texts:
                        if bpy.path.abspath(text.filepath) == self.file_path:
                            with open(self.file_path, 'r', encoding='utf-8') as f:
                                new_text = f.read()
                                text.clear()
                                text.write(new_text)
                                logger(f"Reloaded text block: {text.name}")
                                # try:
                                #     run_script_later(text)
                                # except Exception as e:
                                #     logger("Error while running '" + self.file_path + "':\n", logger(traceback.format_exc()))
                                #     show_warning(str(e))
                            break
                except Exception as e:
                    logger(f"Error reloading {self.file_path}: {e}")

    new_files = [f for f in find_external_python_files() if f not in _watched_files]

    if not new_files:
        logger("No file found...")
        return

    logger("Watching the following new files:\n" + "\n - ".join(new_files))

    for file_path in new_files:
        directory = os.path.dirname(file_path)
        handler = ReloadOnChangeHandler(file_path)
        observer = Observer()
        observer.schedule(handler, path=directory, recursive=False)
        observer.start()
        _observers.append(observer)
        _watched_files.add(file_path)
        logger(f"Now watching: {file_path}")


def watch_for_new_files():
    try:
        logger("Scanning new files...")
        start_file_watchers()
    except Exception as e:
        logger("Error checking for new files:", e)

    return 5.0  # re-run this function every 5 seconds


def _delayed_start():
    try:
        # if bpy.data.filepath:
        logger("Running watcher on currently opened file.")
        start_file_watchers()
    except Exception as e:
        logger("Delayed start failed:", e)
        return 5.0
    return None


@persistent
def on_load(dummy):
    logger("on_load triggered. Starting watchers after file load.")
    bpy.app.timers.register(_delayed_start, first_interval=5.0)


def register():
    logger("Add-on enabled.")
    if ensure_watchdog_installed():
        logger("watchdog loaded with success")

        if on_load not in bpy.app.handlers.load_post:
            bpy.app.handlers.load_post.append(on_load)
            logger("Will load shortly...")

        try:
            logger("Registering _delayed_start")
            bpy.app.timers.register(_delayed_start)
            logger("Registering watch_for_new_files")
            bpy.app.timers.register(watch_for_new_files)
            logger("End of registration")
        except Exception as e:
            logger("Timer registration failed:", e)


def unregister():
    logger("Add-on disabled.")
    if on_load in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.remove(on_load)
        logger("Removed on_load handler")
    for obs in _observers:
        obs.stop()
        obs.join()
    _observers.clear()
    _watched_files.clear()
    logger("Closed all observers")

