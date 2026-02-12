#!/usr/bin/env python3
"""Transform Allen Institute CCFv3 imaging data to standard coordinate spaces.

This script converts neuroimaging volumes between different coordinate systems
and orientations commonly used in mouse brain imaging. It handles conversion
from Allen PIR coordinates to RAS standard imaging coordinates.
"""

# %% Imports

import argparse
import atexit
import numpy as np
import os
import sys
import tempfile
from pyminc.volumes.factory import volumeFromFile, volumeFromDescription
from subprocess import call


def is_minc(infile):
    """Check if file is in MINC format.
    
    Args:
        infile: Path to input file
        
    Returns:
        bool: True if file extension is .mnc
    """
    return os.path.splitext(infile)[1] == ".mnc"


def itk_convert(infile, outfile):
    """Convert between different neuroimaging formats using ITK.
    
    Args:
        infile: Path to input file
        outfile: Path to output file
    """
    call(["itk_convert", "--clobber", infile, outfile])


def remove_temp_files(file_list):
    """Clean up temporary files on exit.
    
    Args:
        file_list: List of temporary file objects to remove
    """
    for f in file_list:
        # Safely obtain the filename associated with this temporary file object.
        filename = getattr(f, "name", None)
        if not filename:
            # Nothing to clean up if there is no filename.
            continue

        # Ensure the file handle is closed, ignoring errors if already closed.
        try:
            if not getattr(f, "closed", False):
                f.close()
        except Exception:
            # Best-effort cleanup; continue to attempt unlinking the path.
            pass

        # Explicitly remove the temporary file from disk if it still exists.
        try:
            if os.path.exists(filename):
                os.remove(filename)
        except FileNotFoundError:
            # File was already removed; nothing more to do.
            pass
def reorient_to_standard(dat):
    """Reorient data from PIR to RAS orientation.
    
    Args:
        dat: Input numpy array
        
    Returns:
        np.ndarray: Reoriented array
    """
    dat = np.rot90(dat, k=1, axes=(0, 2))
    dat = np.rot90(dat, k=1, axes=(0, 1))

    shape = dat.shape
    dat = np.ravel(dat)
    dat = np.reshape(dat, shape)

    return dat


def do_nothing(dat):
    """Identity function for data that doesn't need reorientation.
    
    Args:
        dat: Input numpy array
        
    Returns:
        np.ndarray: Unmodified input array
    """
    return dat


# %% Argument parsing

parser = argparse.ArgumentParser(
    description='Orient Allen Institute CCFv3 imaging data to standard spaces'
)
parser.add_argument('infile', type=str, help='Input volume file path')
parser.add_argument('outfile', type=str, help='Output volume file path')
parser.add_argument('--tmpdir', '-t', default='/tmp', type=str,
                    help='Temporary directory (default: /tmp)')
parser.add_argument('--voxel_orientation', '-v', default='RAS', type=str,
                    choices=['RAS', 'PIR'],
                    help='Voxel orientation: "RAS" (default) or "PIR"')
parser.add_argument('--world_space', '-w', default='MICe', type=str,
                    choices=['MICe', 'CCFv3', 'stereotaxic'],
                    help='World space: "MICe" (default), "CCFv3", or "stereotaxic"')
parser.add_argument('--expansion_factor', '-x', default=1.0, type=float,
                    help='Factor to expand volume by (e.g., 1000 to convert mm to um)')
parser.add_argument('--volume_type', default=None, type=str,
                    help='Volume type (default: from input file)')
parser.add_argument('--data_type', default=None, type=str,
                    help='Data type (default: from input file)')
parser.add_argument('--labels', action='store_true',
                    help='Treat as label volume')
parser.add_argument('--clobber', default=False, action='store_true',
                    help='Overwrite existing output file')

args = parser.parse_args()

# %% Validate inputs

if not os.path.exists(args.infile):
    print(f"Error: Input file not found: {args.infile}", file=sys.stderr)
    sys.exit(1)

if os.path.exists(args.outfile) and not args.clobber:
    print("Error: Output file exists! Use --clobber to overwrite", file=sys.stderr)
    sys.exit(1)

# %% Coordinate definitions

# Centers are listed as x,y,z (in um in CCFv3 coordinates)
centers_RAS = {
    "MICe": [5700, 7900, 2700],
    "CCFv3": [0, 13200, 8000]
}
centers_PIR = {
    "MICe": [5300, 5300, 5700],
    "CCFv3": [0, 0, 0]
}

# Direction cosines
direction_cosines_RAS = {
    "MICe": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "CCFv3": [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]
}

direction_cosines_PIR = {
    "MICe": [[0, -1, 0], [0, 0, -1], [1, 0, 0]],
    "CCFv3": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
}

# Map arguments to functions/dicts
map_voxel_orientations = {
    "RAS": reorient_to_standard,
    "PIR": do_nothing
}

map_centers = {
    "RAS": centers_RAS,
    "PIR": centers_PIR
}

map_dir_cosines = {
    "RAS": direction_cosines_RAS,
    "PIR": direction_cosines_PIR
}

# Map volume sizes to resolutions (in microns)
size_10 = 1320 * 800 * 1140
size_25 = 528 * 320 * 456
size_50 = 264 * 160 * 228
size_100 = 132 * 80 * 114
size_200 = 58 * 41 * 67
map_resolutions = {
    size_10: 10,
    size_25: 25,
    size_50: 50,
    size_100: 100,
    size_200: 200
}

# %% Preprocessing

temporary_files = []

if not is_minc(args.infile):
    tmp = tempfile.NamedTemporaryFile(
        dir=args.tmpdir,
        prefix='mousetools-tmpconv-',
        suffix='.mnc',
        delete=False
    )
    infile = tmp.name
    itk_convert(args.infile, infile)
    temporary_files.append(tmp)
else:
    infile = args.infile

# %% Main script

try:
    # Read volume
    vol = volumeFromFile(infile)

    # Detect resolution
    if vol.data.size not in map_resolutions:
        print(f"Error: Unrecognized volume size: {vol.data.size}", file=sys.stderr)
        print(f"Expected one of: {list(map_resolutions.keys())}", file=sys.stderr)
        sys.exit(1)

    res = map_resolutions[vol.data.size]

    # Voxel orientation
    if args.voxel_orientation in map_voxel_orientations:
        new_data = map_voxel_orientations[args.voxel_orientation](vol.data)
    else:
        print(f"Error: Invalid voxel orientation: {args.voxel_orientation}", file=sys.stderr)
        sys.exit(1)

    # Validate world space
    if args.world_space not in map_centers[args.voxel_orientation]:
        print(f"Error: Invalid world space: {args.world_space}", file=sys.stderr)
        sys.exit(1)

    # World coordinate system
    centers = [
        args.expansion_factor * c / 1000
        for c in map_centers[args.voxel_orientation][args.world_space]
    ]
    steps = [args.expansion_factor * res / 1000] * 3
    xdc = map_dir_cosines[args.voxel_orientation][args.world_space][0]
    ydc = map_dir_cosines[args.voxel_orientation][args.world_space][1]
    zdc = map_dir_cosines[args.voxel_orientation][args.world_space][2]

    # Types
    vtype = vol.volumeType if args.volume_type is None else args.volume_type
    dtype = vol.dtype if args.data_type is None else args.data_type
    # Use True if user provided --labels flag, otherwise use volume's default setting
    labels = True if args.labels else vol.labels

    # %% Data output

    if not is_minc(args.outfile):
        tmp = tempfile.NamedTemporaryFile(
            dir=args.tmpdir,
            prefix='mousetools-tmpconv-',
            suffix='.mnc',
            delete=False
        )
        outfile = tmp.name
        temporary_files.append(tmp)
    else:
        outfile = args.outfile

    outvol = volumeFromDescription(
        outputFilename=outfile,
        dimnames=["zspace", "yspace", "xspace"],
        sizes=new_data.shape,
        starts=[-c for c in reversed(centers)],
        steps=[s for s in reversed(steps)],
        x_dir_cosines=xdc,
        y_dir_cosines=ydc,
        z_dir_cosines=zdc,
        volumeType=vtype,
        dtype=dtype,
        labels=labels
    )
    outvol.data = new_data
    outvol.writeFile()
    outvol.closeVolume()

    # %% Postprocessing

    if not is_minc(args.outfile):
        itk_convert(outfile, args.outfile)

    print(f"Successfully converted {args.infile} to {args.outfile}")

except Exception as e:
    print(f"Error during conversion: {e}", file=sys.stderr)
    sys.exit(1)

# %% On exit

atexit.register(remove_temp_files, file_list=temporary_files)

