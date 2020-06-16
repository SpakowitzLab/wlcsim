.. _tips:

Tips
####

This page contains usefull tips that we have found usefull for using this
codebase.  These may be helpful for your workflow.

Running Remotely
================

When running wlcsim remotely you may want to use tmux so that you don't get
logged off while wlcsim is running.  More documentation on this comming
shortly.

The defines.inc file
====================

You will likely end up with many different defines.inc files.  Generally I
keep a few examples around in `input/exampledefines`.  One of the best tools for comparing defines files is `vim -d defines.inc path/other_defines.inc`.  Look up vim diff online for more shortcuts when using this, they are really helpfull.

Each setting comes with a description.  But you will likely run into cases where the description isn't enough.  Because the settings are named with uinque strings, it is easy to grep the code base for them.  For example:

grep -r "WLC_P__WHATEVERSETTING" src/

Running / storing simulations
=============================

Genearally I coppy the *entire* code base around when I want to run it again with different settings.  The code is small compaired to the output so this donesn't really take up much extra space.  I recompile every time I run with `make clean`.  Since compiling `wlcsim` is fast compaired to running you might as well.

Compiling options
=================

When debugging you likely want to run with `FCFLAGS = ${PEDANTICFLAGS}` in `Makefile`.  When running long simulations you likely want to switch to `FCFLAGS = ${FASTFLAGS}`.  Be sure to `make clean` after changing this setting.

Trouble Shooting
================

When running Monte Carlo many different runtime errors will lead to the integrated energy differing from the energy when it is calculated again from scratch. You'll want to moniter the error output for this warning. This check is one of the best features of this codebase when it comes to code testing, trust us you want it! Here are some tips for if this should occure:

- Turn off some moves and not others.  If the error only occures on some subset that can be informative.

- Are the energies really big?  This may have tripped this error.  You likely don't acctually want them to be that big.

- Which energy type is not matching?  Are there multiple that arn't matching.

If the text in the error message appears to be refering to text that isn't actually there.

- The text may have been replaced by the precomilier to the values givin in `src/defines.inc`.  This probably means there's an error with one of your settings.
