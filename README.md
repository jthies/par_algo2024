# Parallel Algorithms 2024/2025

This repository contains material for the programming tasks of MasterMath course "parallel algorithms" 24/25.
It is a dynamically extended tutorial that will take you from first steps (installing the Fortran compiler and libraries)
to being able to write parallel programs and run them on a supercomputer yourself.

If you have never worked with ``git``, it may be helpful to do this [basic tutorial](https://product.hubspot.com/blog/git-and-github-tutorial-for-beginners).

# Working with the repository

In a terminal (either on your own computer or a cluster like DelftBlue),
first clone the repository to get a local copy:
```bash
git clone https://github.com/jthies/par_algo2024
```

Whenever we announce an update via the [elo course page](https://elo.mastermath.nl/course/view.php?id=1065),
you need to update your local copy:
```bash
git pull
```
If you have local changes to files, this may result in a message:
```bash
error: cannot pull with rebase: You have unstaged changes.
error: Please commit or stash them.
```
The ``git status`` command shows you what files have local changes.
You can run ``git stash`` to move these changes out of the way for the pull.
If you then want to re-do them (after the ``pull``), simply type ``git stash pop``.

**Note:** You are free to commit your work to the local repository or a fork on github,
but you cannot push to the repo itself. Committing your work locally makes sure you have
the history and a backup in case you break your programs and is therefore highly recommended.
Having local commits will typically not cause conflicts with the ``git pull`` as we only
add new files for new tasks.

# 0. Basic programming environment

We recommend that you set up a working environment on your own computer.
During the course, you will also get access to the DelftBlue supercomputer,
where you do not need to install additional software.

The programming language used throughout the course is **Fortran 2018**. It has native support for
parallel programming, so in principle all you will need is a Fortran compiler.
We recommend using the GNU compiler ``gfortran``. To enable actual parallel execution of your code,
it relies on the OpenCoarrays library, which you have to install separately:
See [this page](http://www.opencoarrays.org/) for instructions on how to do this.

On the DelftBlue supercomputer, you do not need to install anything. Instead, you should load the required modules
and you're good to go:
```bash
module load 2024r1 openmpi opencoarrays
```
## Makefiles

Compiling a program (especially if it consists of multiple source files)
can require (a series of) commands that are hard to memorize and specific
to the programming environment. The Linux command ``make`` allows you to put the
rules for "building" your program into a file called ``Makefile``, and then uses
these rules to create the program from the source code.

In this repository you will find a fully functional ``Makefile`` geared towards
use on DelftBlue. On your own system, you need to set the environment variable ``OPENCOARRAYS_ROOT``
to the directory where the library is installed. On my Ubuntu laptop, I installed the library using
```bash
sudo apt install opencoarays-dev libcaf-openmpi-3t64
```
and in the file ``.bashrc`` (which is executed at the start of every terminal session), added the line
```bash
export OPENCOARRAYS_ROOT=/usr/lib/x86_64-linux-gnu/open-coarrays/openmpi/
```
On DelftBlue, it is sufficient to load the module as desribed above.

# 1. Hello World!

You find our first program in ``main_hello.f08``.
You can compile it using:
```bash
make main_hello.x
```
To run it with 4 processes, use:
```bash
mpirun -np 4 ./main_hello.x
```
Look at the program and output it produces.
Do you notice anything unexpected?

## Task: re-order the output

Can you adapt the program so that all "Hello" lines are printed before the "Goodbye" lines?

# 2. Dot Product

The second program is slightly more complicated. It consists of two source files: ``dotprod.f08``
and ``main_dotprod.f08``. The former is a Fortran __module__ that provides several implementations
of a function to compute the inner product of two distributed vectors (coarrays). The latter contains
a program that will run these variants and print timing results.

- Look at both source files and understand (at least roughly) how this program works, and what the different
variants are.
- compile and run the program "as is" using
```bash
make main_dotprod.x
mpirun -np 4 ./main_dotprod.x 1000000
```
You can adapt the number of processes if you have more than 4 cores available. In this example we use one million vector elements in total,
you can change this as well.

## Task: Implement the ``gatherbcast`` variant

You may have noticed that the ``gatherbcast`` variant is not implemented yet.
Do this in ``dotprod.f08`` and run the program again to test and time it.
The way the variant works is sketched in the comments of that function.

## Task: Implement your own variant

Implement at least one additional variant of your own. You can also try to generalize the "butterfly"
variant, which currently only works if the number of processes is a power of 2.

Run the main program for the maximum number of processes available on your system and different values of the vector length N.
Which variant is the best?

