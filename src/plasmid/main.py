#!/usr/bin/env python
# Main function to execute commandline utilities

# system and time
import argparse

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('-v', dest='verbose', action='store_true', help='print verbose info')
    parser.add_argument('-c', dest='config', action='store_true', help='print verbose info')
    subparser = parser.add_subparsers(title='subcommands', dest='subcommand')
    
    # parse arguments
    args = parser.parse_args()
