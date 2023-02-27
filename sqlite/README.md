Setup and use SQLite on Roar:

- Download and uncompress SQLite:

	$ bash download_sql.sh

- Compile and install SQLite:

	$ cd sqlite-autoconf-3410000
	$ ./configure --prefix=/storage/work/<user_id>
	$ make && make install

- Add `/storage/work/<user_id>/bin` to your path (if not already done):

	$ echo "PATH=$PATH:/storage/work/<user_id>/bin" >> ~/.bashrc


References:
 - https://falksangdata.no/wp-content/uploads/2021/04/Opensource.com_cheat_sheet_sqlite_3.pdf
 - https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html
