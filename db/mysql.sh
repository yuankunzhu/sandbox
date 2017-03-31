# load local csv via cmd
mysql -u username -p password --local-infile scrapping -e "LOAD DATA LOCAL INFILE 'CSVname.csv'  INTO TABLE table_name  FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'"

brew services stop postgresql
brew services start postgresql
ps auxwww | grep postgres