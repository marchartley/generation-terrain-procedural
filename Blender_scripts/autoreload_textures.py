import bpy
import os
import time
import sys
import glob
import shutil


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


# Function to collect texture information
def collect_textures():
    texture_mapping = {}
    image_paths = {}
    image_mapping = {}
    for texture in bpy.data.textures:
        if texture.type == 'IMAGE' and texture.image and texture.image.filepath:
            path = bpy.path.abspath(texture.image.filepath)
            texture_mapping[path] = texture.name
            image_mapping[texture.name] = texture.image.name
            if path not in image_paths:
                image_paths[path] = os.stat(path).st_mtime
    return image_paths, texture_mapping, image_mapping


image_paths, texture_mapping, image_mapping = collect_textures()


def applyRender(filename):
    file = f"/media/marc/Data/NN Datasets/1/noise_renders/Render-{filename}"
    bpy.context.scene.render.filepath = file
    if os.path.exists(file + ".png"):
        print(f"File '{file}' already exists")
        print("Continuing")
        return
    scale = 1 / 2
    bpy.context.scene.render.resolution_x = int(1920 * scale)  # = 500 #perhaps set resolution in code
    bpy.context.scene.render.resolution_y = (1080 * scale)  # = 500
    blockPrint()
    bpy.ops.render.render(write_still=True)
    enablePrint()


# Function to check for updates in image files
def check_file_updates(verbose=True):
    changesOccurs = False
    for path, last_modified in image_paths.items():
        try:
            current_modified = os.stat(path).st_mtime
            if current_modified != last_modified:
                image_paths[path] = current_modified  # Update last modified time
                update_texture_in_blender(path, verbose)
                changesOccurs = True
        except FileNotFoundError:
            continue  # File not found, skip this check
    if changesOccurs:
        if verbose:
            print("Some changes!")
        # applyRender()
    return 0.5  # Check every second


# Function to update texture in Blender
def update_texture_in_blender(image_path, verbose=True):
    texture_name = texture_mapping[image_path]
    texture = bpy.data.textures.get(texture_name)
    if texture and texture.type == 'IMAGE' and texture.image:
        texture.image.reload()
        refresh_modifier(texture_name, verbose)
        refresh_materials_using_texture(texture_name, verbose)
        if verbose:
            print(f"Updated texture in Blender: {texture_name}")


# Refresh the displace modifier associated with the texture
def refresh_modifier(texture_name, verbose=True):
    for obj in bpy.data.objects:
        for modifier in obj.modifiers:
            try:
                if modifier.texture and modifier.texture.name == texture_name:
                    # Reassign the texture to refresh the modifier
                    modifier.texture = None
                    modifier.texture = bpy.data.textures[texture_name]
                    if verbose:
                        print(f"Refreshed modifier for {obj.name}")
            except Exception as e:
                pass


# Refresh materials that use the updated texture
def refresh_materials_using_texture(texture_name, verbose=True):
    for mat in bpy.data.materials:
        if mat.use_nodes:
            for node in mat.node_tree.nodes:
                if node.type == 'TEX_IMAGE' and node.image and node.image.name == image_mapping[texture_name]:
                    mat.node_tree.nodes.update()  # Trigger nodes update
                    mat.update_tag()
                    if verbose:
                        print(f"Refreshed material shader nodes for {mat.name}")
                    toggle_viewport_shading()


def toggle_viewport_shading():
    # Loop through all areas in all windows
    for window in bpy.context.window_manager.windows:
        screen = window.screen
        for area in screen.areas:
            if area.type == 'VIEW_3D':  # Check if the area is a 3D viewport
                # Get the active space data for the area
                space = area.spaces.active
                # Save the current shading type
                current_shading = space.shading.type
                # Toggle to wireframe, then back to the original shading to force refresh
                space.shading.type = 'WIREFRAME'
                area.tag_redraw()  # Request redraw of the area
                # Reset the shading type
                space.shading.type = current_shading


def renderSequence():
    toRenderList = sorted(glob.glob("/media/marc/Data/NN Datasets/1/noise_heightmaps/*.png"))
    destination = "/media/marc/Data/NN Datasets/1/result_height.png"
    shutil.copyfile(destination, destination + ".tmp")
    print(toRenderList)
    for i, file in enumerate(toRenderList):
        filename = os.path.basename(file)
        print(f"Rendering ({i + 1}/{len(toRenderList)})... ", end="", flush=True)
        shutil.copyfile(file, destination)
        check_file_updates(verbose=False)
        applyRender(f"{filename}")
        print(f"Done: '{filename}' -> '{bpy.context.scene.render.filepath}'")
    shutil.copyfile(destination + ".tmp", destination)
    check_file_updates()
    print("End.")


if __name__ == "__main__":
    # Register the timer function to check for file updates
    bpy.app.timers.register(check_file_updates)
    # renderSequence()
    # applyRender()


