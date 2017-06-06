import logging
import argparse
import ConfigParser
import sqlalchemy


def main():
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('configfile')
    args = parser.parse_args()
    configs = ConfigParser.ConfigParser()
    configs.read(args.configfile)
    postgres = 'postgresql://{}:{}@{}/{}'.format(
        configs.get('harvestdb', 'username'),
        configs.get('harvestdb', 'password'),
        configs.get('harvestdb', 'hostname'),
        configs.get('harvestdb', 'database'))
    engine = sqlalchemy.create_engine(postgres)
    print engine
    engine.connect()
    # print sqlalchemy.MetaData()
    engine.dispose()

if __name__ == '__main__':
    main()
