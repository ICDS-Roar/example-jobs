# Table of Contents

* [Adding new example jobs](#adding-new-example-jobs)
* [Updating example jobs](#updating-example-jobs)
* [Opening issues](#opening-issues)
* [Suggesting enhancements](#suggesting-enhancements)

# Adding new example jobs

## How to submit an example job

To submit an example job to this repository, simply create a fork of this repository. In that fork, create a branch named after the example job that you wish to submit. Once you complete your example job, open a pull request to main repository. However, before you open the pull request, make sure you follow these steps:

1. Build and test your job.
2. Write a README explaining how to set up your job, how to submit your job, and the output the user should expect to see when the example job completes.
3. Write a bash script named `setup.sh` that handles installing dependencies for the job.

After you open the pull request, we will review the pull request, and we will determine if it should be merged or not! If we deny the merge, we will comment why.

## Criteria that your example job needs to meet

In order to have an example job accepted into this repository, the job must meet the following criteria:

1. The example job's name is indicative of the implementation langauge or application.
2. The example job can successfully run on the open queue.
3. The example job takes under an hour to complete.
4. The example job's environment is reproducible.
5. The example job follows good development practices.

Below, we have expanded upon the criteria that your example job must meet.

### Name is indicative of implementation language or application

As a general rule of thumb, name your job after the implementation langauge/application that you are using. For example, if you are using an example job with Python, Numpy, and Pandas, name the job `python-numpy-pandas`. The same logic applies for jobs that rely on an application. If you have a job that uses Abaqus, name it `abaqus`. Do not get us wrong, we love cool names, however, we need to realize that not everyone might understand what the name you gave to your job means.

### Can successfully run on the open queue

This one is self explanitory, but the reason for this requirement is because every user on the the cluster has access to the open queue. The same cannot be said for private allocations. Therefore, your job should as accessible as possible to the entire user community! Hence the requirement for it to need to successfully run on the open queue.

### Takes under an hour to complete

The point of the jobs here is to demonstrate the general functionality and the ins and outs of submitting batch jobs on Roar. Everyone likes jobs that demonstrate functionality in a clear, concise way in a short amount of time. Therefore, it is neccessary for jobs to take under an hour to complete. To be honest, under 10 minutes is even better!

### The environment is reproducible

Reproducibility is king. For your example job, make sure that you are documenting all the changes that you are making to your environment on Roar. The goal of this is so that another user can easily recreate the environment you used and successfully submit their own copy of the job. If you do not provide instructions for how to setup your job, users will end up getting frustrated and passing your job up!

### Good development practices

Here are some general bullet points that you need to keep in mind when developing your example job:

* Include a `.gitignore` file for the language/application you are using.
  
  * As good practice, it is wise to include `.gitignore` files for the programming language that you are using. C is a popular language in HPC, but you really should not be including the compiled executables in your job (it also makes a mess of the git repository). Leave the compilation up to the user that wants to use your example (just include the compilation instructions you know)!
  * You can consult this GitHub repository for some helpful `.gitignore` files: https://github.com/github/gitignore

* Include a README

  * Everyone loves READMEs! They are super helpful for providing users with the neccessary info they need. In your README, you should include set up instructions, how to submit jobs, and the expected output the user should see when the example job completes.

* Follow a programming language's standards

  * This is just a general practice that you should follow in your coding adventures. If you write your own code to follow a language's standards, it will be easier for others to read your code. Also, the same concept applies to dependency management. You should be installing external packages in a standardized fashion. For example, if you are using a language like Python, you should be installing pip packages by either using a `requirements.txt` or `setup.py` file.

* Automate set up if you can

  * In order to help minimize the time that users are spending on setting up an example job, you should just create a shell script that handles the set up for the user. Typically, on Roar, this is just a bash script named `setup.sh`.

* Verify that your example job is practical

  * We know that sometimes people do not want to ask themselves this question, but it is important to ask yourself if other users of Roar will find you example job helpful. We welcome all contributions, but we need to ensure that the example jobs will be helpful to the Roar user community as whole. A job demonstrating parallelization in R on Roar is super helpful. A job that shows how to print colored text in Golang? Not so much.

# Updating example jobs

If you would like to help us keep our example jobs up to date, you will need to follow these steps:

1. Create a fork of this repository and in that fork create a branch named after the example job you are updating.
2. Build and test the updated example job.
3. Open a pull request to the main repository from your branch.
4. In that pull request, specify what you changed in the existing example job.
5. Wait for your update to be approved or not.

Generally, we will always approve working updates to example jobs. If your pull request is denied, it might be because we either found an issue while testing your updates, or we found that your update did not meet one or more of the aforementioned criterium. Either way, we will be sure to explain why.

# Opening issues

### Bugs

If you encounter any bugs or any *oddities* while working with one of the example jobs stored in this repository, please open an issue on this repository. In that issue, please include the following sections:

1. Which example job are you trying to submit?
2. The stacktrace of the error you are receiving.

The more information the better. We cannot fix the problem if we do not know how it is being caused. Also, when you open the issue, please label the issue as a **bug**.

### Documentation

If you have any questions about the documentation or anything mentioned there within, please open an issue marked as **documentation**. The goal of the documentation here is that it will be useful in getting you up to speed on a topic. If something is difficult to understand or something is outdated, please let us know! We will be sure to fix/update the relevant documentation.

# Suggesting enhancements

If there is a new example job you would like to see added to this repository or any modifications existing example jobs that you would like to see, please open an issue on this repository. While we cannot promise that every enhancemnet request will be followed through on, we will at least give it a look! Also, when requesting a feature as an issue, please label the issue as a **feature request** or **enhancement**.
