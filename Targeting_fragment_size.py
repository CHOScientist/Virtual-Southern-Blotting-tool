#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integrated Script for Enzyme Distance Calculations and Finalizing Fragment Lengths

Description:
------------
This script merges the functionality of three separate code blocks used for:
1) Reading enzyme positions and computing distances to integration sites.
2) Restructuring or preparing the distance output.
3) Incorporating fragment lengths and directions to produce final distances.

Usage:
------
    1) Calculate distances from enzyme sites:
       python integrated_script.py calc-distances \
           --enzyme-file PICRH_enzyme_site.csv \
           --integration-file integration_site.csv \
           --output closest_enzyme_distances_new.csv

    2) Finalize distances by merging fragment lengths and directions:
       python integrated_script.py finalize-distances \
           --lengths-file fragment_lengths.csv \
           --directions-file directions.csv \
           --distances-file restructured_closest_enzyme_distances_new.csv \
           --output final_distances_with_lengths.csv

Dependencies:
-------------
- Python 3.x
- csv (built-in)
- Any other libraries if needed

Author:
-------
Your Name (your.email@example.com)
"""

import csv
import argparse
import sys


# ============================================================================
# Section 1: Reading and Calculating Distances from Enzyme Positions
# ============================================================================

def read_enzyme_positions(filename):
    """
    Reads enzyme site positions from a CSV file. Each row must have:
    [enzyme_name, chromosome, position].

    Parameters
    ----------
    filename : str
        Path to the CSV file containing enzyme positions.

    Returns
    -------
    dict
        Nested dictionary of the form:
        {
            enzyme_name: {
                chromosome: [sorted positions (int)]
            },
            ...
        }
    """
    enzyme_positions = {}
    with open(filename, 'r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) == 3:
                enzyme, chromosome, position = row
                if enzyme not in enzyme_positions:
                    enzyme_positions[enzyme] = {}
                if chromosome not in enzyme_positions[enzyme]:
                    enzyme_positions[enzyme][chromosome] = []
                enzyme_positions[enzyme][chromosome].append(int(position))
            else:
                print(f"[WARNING] Skipping malformed enzyme row: {row}")

    # Sort positions for each chromosome
    for enzyme in enzyme_positions:
        for chrom in enzyme_positions[enzyme]:
            enzyme_positions[enzyme][chrom].sort()

    return enzyme_positions


def find_closest(positions, target):
    """
    Finds the closest upstream and downstream positions relative to a target.

    Parameters
    ----------
    positions : list of int
        Sorted list of positions for a given chromosome.
    target : int
        The reference position for which to find closest upstream/downstream sites.

    Returns
    -------
    tuple
        (upstream, downstream), where each is an integer or None.
    """
    upstream = None
    downstream = None
    for pos in positions:
        if pos < target:
            upstream = pos
        elif pos > target and downstream is None:
            downstream = pos
            break
    return upstream, downstream


def calculate_distances(enzyme_positions, integration_sites):
    """
    Calculates distances from each integration site to the closest upstream
    and downstream enzyme sites.

    Parameters
    ----------
    enzyme_positions : dict
        Nested dict of enzyme positions created by read_enzyme_positions().
    integration_sites : str
        Path to the CSV file listing integration sites [site_id, chromosome, position].

    Returns
    -------
    list of list
        Each sublist contains [site_id, chromosome, position, enzyme, upstream_dist, downstream_dist, ...].
    """
    results = []
    with open(integration_sites, 'r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) == 3:
                site_id, chromosome, position = row
                position = int(position)
                row_data = [site_id, chromosome, position]

                # For each enzyme, find distances
                for enzyme in enzyme_positions:
                    if chromosome in enzyme_positions[enzyme]:
                        upstream, downstream = find_closest(enzyme_positions[enzyme][chromosome], position)
                        upstream_distance = position - upstream if upstream is not None else "N/A"
                        downstream_distance = downstream - position if downstream is not None else "N/A"
                    else:
                        upstream_distance = downstream_distance = "N/A"

                    # Append enzyme + distances
                    row_data.extend([enzyme, upstream_distance, downstream_distance])

                results.append(row_data)
            else:
                print(f"[WARNING] Skipping malformed integration row: {row}")

    return results


def write_results(results, filename):
    """
    Writes the distance calculation results to a CSV file.

    Parameters
    ----------
    results : list of list
        List of rows generated by calculate_distances().
    filename : str
        Output CSV file path.
    """
    # Collect all unique enzymes (assuming they appear at fixed indices)
    # Example row: [site_id, chrom, pos, EnzymeA, upDistA, downDistA, EnzymeB, ...]
    # The first 3 columns are always site_id, chromosome, position.
    # Then for each enzyme, we have 3 columns: [enzyme_name, upDist, downDist].
    if not results:
        print("[INFO] No results to write.")
        return

    # Identify unique enzyme names by scanning one row
    # The pattern in each row after the first three columns is (enzyme, upstream_dist, downstream_dist)
    row_example = results[0]
    # first 3 = [site_id, chromosome, position]
    # after that = repeated sets of 3
    enzyme_names = []
    chunked = row_example[3:]
    for i in range(0, len(chunked), 3):
        enzyme_names.append(chunked[i])

    # Prepare CSV header
    header = ['IntegrationSite#', 'Chromosome', 'Position']
    for enzyme in enzyme_names:
        header.append(f"{enzyme} UpstreamDist")
        header.append(f"{enzyme} DownstreamDist")

    with open(filename, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for row_data in results:
            writer.writerow(row_data)

    print(f"[INFO] Distance results have been written to {filename}")


# ============================================================================
# Section 2: Finalizing Distances with Fragment Lengths and Directions
# ============================================================================

def read_fragment_lengths(filename):
    """
    Reads fragment lengths from a CSV file with columns:
    [Name, L, H].

    Parameters
    ----------
    filename : str
        Path to the CSV file.

    Returns
    -------
    dict
        {
            'EnzymeName': {'L': int, 'H': int},
            ...
        }
    """
    lengths = {}
    with open(filename, 'r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            # Example row: Name=EnzymeA, L=100, H=200
            enzyme_name = row['Name']
            lengths[enzyme_name] = {
                'L': int(row['L']),
                'H': int(row['H'])
            }
    return lengths


def read_directions(filename):
    """
    Reads direction information from a CSV file with columns:
    [IntegrationSite#, <enzyme>_L=up/down, <enzyme>_H=up/down, etc.].

    Parameters
    ----------
    filename : str
        Path to the CSV file containing direction info.

    Returns
    -------
    dict
        {
            'SiteID': {
                '<enzyme>_L': 'up' or 'down',
                '<enzyme>_H': 'up' or 'down',
                ...
            },
            ...
        }
    """
    directions = {}
    with open(filename, 'r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            # row['IntegrationSite#'], row['EnzymeA_L'] => up/down
            site_id = row['IntegrationSite#']
            directions[site_id] = dict(row)
    return directions


def calculate_final_distances(distances_file, directions, lengths):
    """
    Integrates distance information with fragment lengths and directions
    to compute final distances.

    Parameters
    ----------
    distances_file : str
        CSV file containing restructured distance info
        (e.g., restructured_closest_enzyme_distances_new.csv).
    directions : dict
        Dict of site directions from read_directions().
    lengths : dict
        Dict of fragment lengths from read_fragment_lengths().

    Returns
    -------
    tuple (list, list)
        - A list of dictionaries with final distance calculations.
        - A list of fieldnames for writing to CSV.
    """
    final_distances = []
    dynamic_fieldnames = set(['IntegrationSite#', 'Chromosome', 'Position'])

    with open(distances_file, 'r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            site_id = row['IntegrationSite#']
            new_row = {
                'IntegrationSite#': site_id,
                'Chromosome': row['Chromosome'],
                'Position': row['Position']
            }

            # For each known enzyme from the lengths dictionary
            for enzyme in lengths.keys():
                # We have L and H in the lengths
                for direction_type in ['L', 'H']:
                    direction_key = f"{enzyme}_{direction_type}"
                    # Build the field name for distances: upstream or downstream
                    # Look up in the directions if this site is 'up' or 'down'
                    # e.g., if directions[site_id]['EnzymeA_L'] == 'up'
                    if site_id in directions:
                        site_directions = directions[site_id]
                    else:
                        site_directions = {}

                    # Identify if L/H is associated with upstream or downstream
                    user_direction = site_directions.get(direction_key, '')
                    if user_direction == 'up':
                        dist_key = f"{enzyme} UpstreamDist"
                    elif user_direction == 'down':
                        dist_key = f"{enzyme} DownstreamDist"
                    else:
                        # If not specified or missing, treat as N/A
                        dist_key = None

                    # If dist_key is valid and present in the CSV row
                    if dist_key and dist_key in row and row[dist_key] != 'N/A':
                        base_distance = int(row[dist_key])
                        final_distance = base_distance + lengths[enzyme][direction_type]
                        new_row[direction_key] = final_distance
                    else:
                        new_row[direction_key] = 'N/A'

                    dynamic_fieldnames.add(direction_key)

            final_distances.append(new_row)

    return final_distances, list(dynamic_fieldnames)


def write_final_distances(final_distances, fieldnames, output_filename):
    """
    Writes the final integrated distances to a CSV file.

    Parameters
    ----------
    final_distances : list of dict
        List of dictionaries with the final distance data.
    fieldnames : list of str
        Column names for the output CSV.
    output_filename : str
        Path to write the final CSV output.
    """
    with open(output_filename, 'w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for row in final_distances:
            writer.writerow(row)

    print(f"[INFO] Final distances with lengths have been written to {output_filename}")


# ============================================================================
# Main Entry Point with Subcommands
# ============================================================================

def main():
    """
    Main entry point. Uses argparse to differentiate between:
    1) 'calc-distances'  - For enzyme distance calculations.
    2) 'finalize-distances' - For merging fragment lengths/directions with distances.
    """
    parser = argparse.ArgumentParser(
        description="Integrated script for enzyme distance calculations "
                    "and final distance finalization."
    )

    subparsers = parser.add_subparsers(dest='command', required=True)

    # Subcommand: calc-distances
    calc_parser = subparsers.add_parser(
        'calc-distances',
        help="Calculate distances from enzyme positions to integration sites."
    )
    calc_parser.add_argument(
        '--enzyme-file',
        required=True,
        help="CSV file containing enzyme positions."
    )
    calc_parser.add_argument(
        '--integration-file',
        required=True,
        help="CSV file containing integration sites."
    )
    calc_parser.add_argument(
        '--output',
        default='closest_enzyme_distances_new.csv',
        help="Output CSV file for writing distance results."
    )

    # Subcommand: finalize-distances
    finalize_parser = subparsers.add_parser(
        'finalize-distances',
        help="Finalize distances by merging fragment lengths and directions."
    )
    finalize_parser.add_argument(
        '--lengths-file',
        required=True,
        help="CSV file containing fragment lengths (Name, L, H)."
    )
    finalize_parser.add_argument(
        '--directions-file',
        required=True,
        help="CSV file containing direction info (IntegrationSite#, Enzyme_L=up/down, etc.)."
    )
    finalize_parser.add_argument(
        '--distances-file',
        required=True,
        help="CSV file containing restructured distance info."
    )
    finalize_parser.add_argument(
        '--output',
        default='final_distances_with_lengths.csv',
        help="Output CSV file for final distances."
    )

    args = parser.parse_args()

    # Handle subcommands
    if args.command == 'calc-distances':
        # Calculate distances from enzyme-file to integration-file
        enzyme_positions = read_enzyme_positions(args.enzyme_file)
        results = calculate_distances(enzyme_positions, args.integration_file)
        write_results(results, args.output)

    elif args.command == 'finalize-distances':
        # Merge distances with fragment lengths and directions
        lengths = read_fragment_lengths(args.lengths_file)
        directions = read_directions(args.directions_file)
        final_distances, fieldnames = calculate_final_distances(
            args.distances_file, directions, lengths
        )
        write_final_distances(final_distances, fieldnames, args.output)


if __name__ == "__main__":
    main()
