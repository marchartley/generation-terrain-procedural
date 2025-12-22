import glob
import re

import PIL.Image
import bpy
import json
import os
import random
import mathutils
import shutil
import math
from PIL import Image
import numpy as np
import PIL.ImageFilter
from mathutils import Euler

timestamp = "latest"
timestamp = "2025-06-20__10-32-20"
timestamp = -1

display_objects = False
display_flow = False # Test
display_material = None # "sand"  # None

use_objects_scores = True

recompute_heights = True

use_127_as_ground = False

minimum_flow_force = 0.05

objects_initial_terrain_dimensions = np.array([100, 100])
heightmap_dimensions = np.array([500, 500])

ratio_obj_to_heightmap = heightmap_dimensions / objects_initial_terrain_dimensions

folder = "/home/marc/generation-terrain-procedural/EnvObjRendering/"  # + timestamp + "/"

if timestamp == "latest":
    timestamp = sorted(os.listdir(folder))[-1]
elif isinstance(timestamp, int):
    timestamp = sorted(os.listdir(folder))[timestamp]

folder = folder + timestamp + "/"


# Path to your JSON file
json_file_path = folder + "terrain_saved.json"
heightmap_path = folder + "heightmap_subsided.png"
flow_path = folder + "flowfield_total.png"

pattern = re.compile(r"obj_score_(.*)\.png")

objects_scores_images = {
    pattern.search(os.path.basename(img_path)).group(1): np.array(Image.open(img_path))
    for img_path in glob.glob(folder + "obj_score_*.png")
    if pattern.search(os.path.basename(img_path))
}

def smooth_heightmap(original: Image):
    return original.resize((heightmap_dimensions[0], heightmap_dimensions[1]), resample=PIL.Image.Resampling.BILINEAR).filter(PIL.ImageFilter.GaussianBlur(radius=5.0))

def compute_new_heightmap(resulting_filename):
    height_details_path = folder + "heightmap_surface-constraint.png"
    # Load images as grayscale
    img1 = Image.open(heightmap_path).convert("I;16")
    img2 = Image.open(height_details_path).convert("I;16")

    # Convert to numpy arrays
    base = np.array(img1).astype(np.float32) * (255.0 / 65535.0) + 127
    disp_raw = np.array(img2).astype(np.float32) * (255.0 / 65535.0)
    displacement = (disp_raw - 127) * 2.0
    result = base + displacement
    # Convert back to image and save
    if not use_127_as_ground:
        result = (result - 127) * 1
    clipped = np.clip(result, 0, 255)
    img3 = smooth_heightmap(Image.fromarray(clipped.astype(np.uint8)))
    img3.save(resulting_filename)

    return np.array(img3)


terrain_object = bpy.data.objects["Terrain with EnvObjs"]

if recompute_heights:
    heightmap = compute_new_heightmap("/media/marc/Data/NN Datasets/1/result_height.png")
else:
    shutil.copy(heightmap_path, "/media/marc/Data/NN Datasets/1/result_height.png")
    heightmap = np.array(smooth_heightmap(Image.open("/media/marc/Data/NN Datasets/1/result_height.png").convert("L"))) - 127.0

if display_material:
    material_path = folder + "material_" + display_material.lower() + ".png"
    shutil.copy(material_path, "/media/marc/Data/NN Datasets/1/material_to_display.png")
    nodes = terrain_object.active_material.node_tree.nodes
    for node in nodes:
        if node.name.lower() == "usedistrib":
            node.outputs[0].default_value = 1.0
else:
    nodes = terrain_object.active_material.node_tree.nodes
    for node in nodes:
        if node.name.lower() == "usedistrib":
            node.outputs[0].default_value = 0

custom_property_name = "created_by_script"
health_property_name = "health"

def getHeight(height_map: np.ndarray, x, y):
    _x, _y = x * ratio_obj_to_heightmap[0], y * ratio_obj_to_heightmap[1]
    terrain_offset = height_map[int(_y), int(_x)]# + 127 #(127.0 if recompute_heights else 0)
    return terrain_offset

def cleanup():
    # Remove all objects with the custom property
    for obj in bpy.data.objects:
        if obj.get(custom_property_name):
            bpy.data.objects.remove(obj, do_unlink=True)

    for texture in bpy.data.textures:
        if texture and texture.type == 'IMAGE' and texture.image:
            texture.image.reload()
            # print(f"Reloaded tex: {texture}")

    for mat in bpy.data.materials:
        if mat.use_nodes:
            for node in mat.node_tree.nodes:
                if node.type == 'TEX_IMAGE' and node.image:  # and node.image.name.lower() == "heightmap":
                    # print(f"Reloaded mat: {mat}")
                    mat.node_tree.nodes.update()  # Trigger nodes update
                    mat.update_tag()


def add_objects():
    # Load JSON data
    with open(json_file_path, 'r') as f:
        data = json.load(f)
        if "instances" not in data:
            data["instances"] = []

    # Enable the OBJ import/export addon
    bpy.ops.preferences.addon_enable(module='io_scene_obj')

    # Supported file extensions and their corresponding import operators
    file_extensions = {
        # '.fbx': bpy.ops.import_scene.fbx,
        # '.stl': bpy.ops.import_mesh.stl,
        # '.obj': bpy.ops.import_scene.obj,
        # '.ply': bpy.ops.import_mesh.ply,
        '.gltf': bpy.ops.import_scene.gltf,
        '.glb': bpy.ops.import_scene.gltf
    }

    def normalize_object(obj, asset_name, flip=False):
        if obj.type != 'MESH':
            print(f"Skipping normalization: {obj.name} is not a mesh")
            return None

        # Force evaluation of the mesh
        depsgraph = bpy.context.evaluated_depsgraph_get()
        obj_eval = obj.evaluated_get(depsgraph)
        mesh = obj_eval.to_mesh()

        if not mesh.vertices:
            print(f"Skipping normalization: {obj.name} has no vertices")
            obj_eval.to_mesh_clear()
            return None

        # Calculate bounding box dimensions
        bbox = [v.co for v in mesh.vertices]
        min_corner = [min(axis) for axis in zip(*bbox)]
        max_corner = [max(axis) for axis in zip(*bbox)]
        dimensions = [max_c - min_c for min_c, max_c in zip(min_corner, max_corner)]

        max_dimension = max(dimensions)
        max_dimension = sum(dimensions) / 3.0
        if max_dimension == 0:
            obj_eval.to_mesh_clear()
            return None

        scale_factor = 1.0 / max_dimension
        # obj.scale = [s * scale_factor for s in obj.scale]
        obj.scale = [scale_factor for s in obj.scale]
        if flip:
            obj.scale[2] *= -1.0

        if asset_name.lower().startswith("arch"):
            rot_quat = mathutils.Vector([0, 0, 90]).to_track_quat('-Z', 'Y')
            # arrow.rotation_euler = rot_quat.to_euler()
            obj.rotation_quaternion = rot_quat

        # Clean up
        obj_eval.to_mesh_clear()
        bpy.context.view_layer.update()
        return obj

    # Create a root object to parent all imported objects
    try:
        root_object = bpy.data.objects["Root"]
    except KeyError:
        root_object = bpy.data.objects.new("Root", None)
        bpy.context.collection.objects.link(root_object)

    # Dictionary to store imported assets
    imported_assets = {}

    def import_and_cache_assets(asset_name, asset_folder, collection_name="ImportedAssets"):
        collection = bpy.data.collections.get(collection_name)
        if not collection:
            collection = bpy.data.collections.new(collection_name)
            bpy.context.scene.collection.children.link(collection)

        cached_objects = []

        # Check if this asset already exists
        if any(obj.name.startswith(asset_name + "_") for obj in collection.objects):
            # Reuse existing cached versions
            cached_objects = [obj for obj in collection.objects if obj.name.startswith(asset_name + "_")]
            for o in cached_objects:
                bpy.data.objects.remove(o, do_unlink=True)
            # return cached_objects

        # Import new files and cache them
        supported_files = [f for f in os.listdir(asset_folder) if os.path.splitext(f)[1].lower() in file_extensions]
        for file in supported_files:
            file_path = os.path.join(asset_folder, file)
            imported_object = import_file(file_path, asset_name)
            if imported_object is None:
                print(f"WARNING: '{file_path}' is an empty object!")
                continue
            imported_object.name = f"{asset_name}_{file}"  # Unique name
            collection.objects.link(imported_object)
            bpy.context.collection.objects.unlink(imported_object)  # Remove from default scene collection
            cached_objects.append(imported_object)
            print(f"Imported: {file_path}")

        return cached_objects

    # Function to import a file
    def import_file(file_path, asset_name):
        file_ext = os.path.splitext(file_path)[1].lower()
        file_import_operator = file_extensions[file_ext]

        # Track objects before import
        existing_objects = set(bpy.data.objects)

        # Perform import
        print(f"Importing: {file_path}")
        result = file_import_operator(filepath=file_path)
        print(result)

        # Track new objects
        new_objects = set(bpy.data.objects) - existing_objects
        if not new_objects:
            print(f"WARNING: No objects imported from {file_path}")
            return None

        # Pick the first newly imported object
        imported_object = new_objects.pop()
        print(f"Imported object: {imported_object.name}")

        # Normalize
        if normalize_object(imported_object, asset_name, flip=file_path.endswith("fbx")) is None:
            print(f"WARNING: Imported object from {file_path} is empty after normalization.")
            return None

        imported_object[custom_property_name] = True
        imported_object[health_property_name] = 1.0

        addMixingInShader(imported_object, (0, 0, 0))

        return imported_object

    for asset_name, asset_folder in data["assets"].items():
        if asset_name in data["instances"]:
            print(f"Processing asset: {asset_name}")
            imported_assets[asset_name] = import_and_cache_assets(asset_name, asset_folder)

    # Iterate over each asset and its instances
    min_z, max_z = 1000000, -10000000
    for asset_name, asset_folder in data["assets"].items():
        if asset_name in data["instances"]:
            # Create instances of the imported asset
            for instance in data["instances"][asset_name]:
                if len(imported_assets[asset_name]) == 0:
                    #print(f"WARNING: asset '{asset_name}' has no models")
                    continue
                # Randomly select an imported object for this instance
                original_object = random.choice(imported_assets[asset_name])
                new_instance = original_object.copy()
                new_instance.data = original_object.data.copy()
                bpy.context.collection.objects.link(new_instance)

                # Set transformation
                new_instance.location = instance["position"]
                new_instance.location[2] = (getHeight(heightmap, instance["position"][0], instance["position"][1]) - 127) * 0.095
                min_z, max_z = min(min_z, new_instance.location[2]), max(max_z, new_instance.location[2])
                # new_instance.location[2] = heightmap[int(instance["position"][1]), int(instance["position"][0])] * 0.095
                new_instance.location -= mathutils.Vector([50, 50, 0])
                new_instance.location = [l / s for l, s in zip(new_instance.location, terrain_object.scale)]
                new_instance.location *= 1.2  # Magic number!
                # new_instance.location -= mathutils.Vector([terrain_object.scale[0], terrain_object.scale[1], 0])
                new_instance.location = [new_instance.location[0], -new_instance.location[1],
                                         new_instance.location[2] - 0.1]

                # new_instance.location -= terrain_object.location
                rot = [new_instance.rotation_euler[i] + instance["rotation"][i] for i in range(3)]
                rot[0] += 90
                rot[2] += ((random.random() - 0.5) * 2.0) * 20
                new_instance.rotation_euler = rot
                # new_instance.rotation_euler[2] += ((random.random() - 0.5) * 2.0) * 90
                # rot_quat = new_instance.rotation_euler.to_quaternion()
                # arrow.rotation_euler = rot_quat.to_euler()
                # new_instance.rotation_quaternion = rot_quat
                # new_instance.rotation_euler[0] += 30
                new_instance.scale *= 0.1
                new_instance.scale = [s * sf for s, sf in zip(new_instance.scale, instance["scale"])]

                obj_score = 1.0
                if use_objects_scores:
                    obj_score = (instance["score"] if "score"in instance else objects_scores_images[asset_name][int(instance["position"][1]), int(instance["position"][0])] / 65535.0 if asset_name in objects_scores_images else 1.0)
                    new_instance.scale *= obj_score
                new_instance[health_property_name] = obj_score

                # Parent the instance to the root object
                new_instance.parent = root_object

    # print("Minimum/Maximum height found: ", min_z, " -- ", max_z)
    # Select the root object and make it active
    bpy.context.view_layer.objects.active = root_object
    bpy.ops.object.select_all(action='DESELECT')
    root_object.select_set(True)


def addMixingInShader(obj, color):
    mat = obj.active_material
    if not mat or not mat.use_nodes:
        print("Material or nodes not found.")
    else:
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links

        tex_node = nodes.get("baseColorTexture.tex_img")
        bsdf_node = None

        for node in nodes:
            if node.type == 'BSDF_PRINCIPLED':
                bsdf_node = node
                break

        if tex_node and bsdf_node:
            # Create MixRGB node
            mix_node = nodes.new(type="ShaderNodeMixRGB")
            mix_node.label = "BASE COLOR MIX"
            mix_node.blend_type = 'MIX'
            mix_node.inputs['Color2'].default_value = (*color, 1.0)

            mix_node.location = ((tex_node.location.x + bsdf_node.location.x) / 2,
                                 (tex_node.location.y + bsdf_node.location.y) / 2)

            # Remove any old link
            for link in links:
                if link.from_node == tex_node and link.to_node == bsdf_node and link.to_socket == bsdf_node.inputs[
                    'Base Color']:
                    links.remove(link)
                    break

            # Reconnect links
            links.new(tex_node.outputs['Color'], mix_node.inputs['Color1'])
            links.new(mix_node.outputs['Color'], bsdf_node.inputs['Base Color'])

            # Set up driver on Fac input
            fac_input = mix_node.inputs['Fac']
            if fac_input.is_linked:
                print("Fac input is already linked. Can't add driver.")
            else:
                # Create the Attribute node
                attr_node = nodes.new(type="ShaderNodeAttribute")
                attr_node.name = health_property_name
                attr_node.label = health_property_name
                attr_node.attribute_name = health_property_name
                attr_node.location = (mix_node.location.x - 200, mix_node.location.y - 100)

                # Connect Attribute 'Fac' output to MixRGB 'Fac' input
                links.new(attr_node.outputs["Fac"], mix_node.inputs["Fac"])

                print("Attribute node connected to MixRGB.Fac.")
        else:
            print("Required nodes not found.")


def add_flow(image_name):
    # Parameters
    # image_name = "/home/marc/generation-terrain-procedural/EnvObjRendering/2025-06-19__23-26-01/flowfield_total.png"
    spacing = 3  # skip every 4 pixels
    arrow_length = spacing  # 1.5
    arrow_radius = 1.5

    # Load image (make sure it's packed or in file path)
    image = bpy.data.images.get(image_name)
    if not image:
        image = bpy.data.images.load(image_name)
    image.colorspace_settings.name = 'Non-Color'

    width = image.size[0]
    height = image.size[1]
    pixels = image.pixels[:]  # flat list: [R, G, B, A, R, G, B, A, ...]

    for mat in bpy.data.materials:
        if mat.use_nodes:
            for node in mat.node_tree.nodes:
                if node.type == 'TEX_IMAGE' and node.image and node.image.name == "Flow":
                    node.image = image
                    mat.node_tree.nodes.update()  # Trigger nodes update
                    mat.update_tag()

    collection_name = "Arrow_flows"
    collection = bpy.data.collections.get(collection_name)
    if not collection:
        collection = bpy.data.collections.new(collection_name)
        bpy.context.scene.collection.children.link(collection)

    # Get/create base arrow
    base_arrow = bpy.data.objects.get("BaseArrow")
    if not base_arrow:
        # Create shared material with Object Info node for per-object color
        mat = bpy.data.materials.get("ArrowMaterial")
        if not mat:
            mat = bpy.data.materials.new(name="ArrowMaterial")
            mat.use_nodes = True

            nodes = mat.node_tree.nodes
            links = mat.node_tree.links
            nodes.clear()

            # Nodes
            output = nodes.new(type="ShaderNodeOutputMaterial")
            bsdf = nodes.new(type="ShaderNodeBsdfPrincipled")
            object_info = nodes.new(type="ShaderNodeObjectInfo")

            # Link object color to base color
            links.new(object_info.outputs['Color'], bsdf.inputs['Base Color'])
            links.new(bsdf.outputs['BSDF'], output.inputs['Surface'])

        shaft_length = 1.2
        head_length = 0.3
        shaft_radius = 0.1
        head_radius = 0.2
        angle_deg = 45
        angle_rad = math.radians(angle_deg)
        direction = mathutils.Vector((math.cos(angle_rad), math.sin(angle_rad), 0))
        # Create base arrow components
        # Shaft (cylinder)
        bpy.ops.mesh.primitive_cylinder_add(vertices=12, radius=shaft_radius, depth=shaft_length)
        shaft = bpy.context.object
        shaft.name = "BaseArrowShaft"
        shaft.location = (0, 0, shaft_length / 2)

        # Head (cone)
        bpy.ops.mesh.primitive_cone_add(vertices=12, radius1=head_radius, depth=head_length)
        head = bpy.context.object
        head.name = "BaseArrowHead"
        head.location = (0, 0, shaft_length + head_length / 2)

        # Join shaft and head into one object
        bpy.ops.object.select_all(action='DESELECT')
        shaft.select_set(True)
        head.select_set(True)
        bpy.context.view_layer.objects.active = shaft
        bpy.ops.object.join()
        base_arrow = shaft
        base_arrow.name = "BaseArrow"

        base_arrow.hide_viewport = True
        base_arrow.hide_render = True
        base_arrow.cycles_visibility.shadow = False

        # Assign material
        base_arrow.data.materials.append(mat)

    # Remove old arrows
    for obj in bpy.data.objects:
        if obj.name.startswith("Arrow_"):
            bpy.data.objects.remove(obj, do_unlink=True)

    # Helper to decode flow from RGB
    def rgb_to_vector(r, g, b):
        def map_channel(val):
            return (val - 0.5) * 2  # map [0,1] to [-1,1]

        return mathutils.Vector((
            -map_channel(r),
            map_channel(g),
            map_channel(b)
        ))  # .normalized()

    terrain_image = bpy.data.images.get("heightmap")

    terrain_width = terrain_image.size[0]
    terrain_height = terrain_image.size[1]
    terrain_pixels = terrain_image.pixels[:]  # flat list: [R, G, B, A, R, G, B, A, ...]

    # Create arrows from image
    for y in range(0, height, spacing):
        for x in range(0, width, spacing):
            index = (y * width + x) * 4  # 4 channels per pixel
            r = pixels[index]
            g = pixels[index + 1]
            b = pixels[index + 2]

            r, g, b = r, g, b

            dx, dy = spacing, spacing
            i00 = (min(max(int(y - dy), 0), height - 1) * width + min(max(int(x - dx), 0), width - 1)) * 4
            i10 = (min(max(int(y + dy), 0), height - 1) * width + min(max(int(x - dx), 0), width - 1)) * 4
            i11 = (min(max(int(y + dy), 0), height - 1) * width + min(max(int(x + dx), 0), width - 1)) * 4
            i01 = (min(max(int(y - dy), 0), height - 1) * width + min(max(int(x + dx), 0), width - 1)) * 4
            """
            start_height = terrain_pixels[i0]
            end_height = terrain_pixels[i1]
            b = 0.5 + (end_height - start_height)
            print(start_height, end_height, b)"""

            # Convert to direction vector
            direction = rgb_to_vector(r, g, b)
            length = direction.length
            if length < minimum_flow_force:
                continue  # skip zero-flow

            # Duplicate base arrow
            arrow = base_arrow.copy()
            arrow.data = base_arrow.data.copy()
            arrow.name = f"Arrow_{x}_{y}"
            bpy.context.collection.objects.link(arrow)

            terrain_offset = max([terrain_pixels[indx] for indx in [i00, i10, i01, i11, index]])  # / 5.0

            # Position and rotation
            pos = mathutils.Vector((x, y, 0))
            arrow.location = pos + direction * length * (arrow_length / 2)
            arrow.location -= mathutils.Vector([50, 50, 0])
            arrow.location[2] = terrain_offset * 20 - 20
            arrow.location = [l / s for l, s in zip(arrow.location, terrain_object.scale)]
            arrow.location[2] += 3.5

            arrow.scale = [2.0 / s for s in terrain_object.scale]
            arrow.scale[2] *= length ** 0.1

            # Align Z-axis of arrow to flow direction
            rot_quat = direction.to_track_quat('-Z', 'Y')
            arrow.rotation_euler = rot_quat.to_euler()

            # Optional: color per arrow
            arrow.color = (r, g, b, 1)

            arrow.hide_viewport = False
            arrow.hide_render = False

            arrow[custom_property_name] = True

            # collection.objects.link(arrow)
            # bpy.context.collection.objects.unlink(arrow)  # Remove from default scene collection


cleanup()

if display_objects:
    add_objects()

if display_flow:
    add_flow(flow_path)

print(f"All ready for '{folder}'!")