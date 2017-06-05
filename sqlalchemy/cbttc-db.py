import logging
import argparse
import ConfigParser


def main():
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('configfile')
    args = parser.parse_args()
    configs = ConfigParser.ConfigParser(args.configfile)
    configs.read(configs)

if __name__ == '__main__':
    main()
