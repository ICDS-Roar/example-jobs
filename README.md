# Pending Brnach

This branch contains broken or archival examples which need testing or other modifications before they are ready to be made live.

To migrate examples to the main branch:

Check out the main branch

	$ git checkout main

Create a new branch for your edits, using a meaningful name such as one that includes the software title and date

	$ git checkout -b software-edits-19790101

Checkout any relevent files on the `pending` branch

	$ git checkout pending -- path/to/file

Make necessary edits and test your changes to ensure they run as expected on Roar. Then make a commit with a message that details the changes made including the software title(s)

	$ git add .
	$ git commit -m "message describing the changes"

Finally, make pull request on the GitHub repository website
