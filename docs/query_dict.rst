.. highlight:: catalog

**********
Query dict
**********

This document describes how to generate a query
*dict* for use with a `specdb` database.

Basics
======

Aside from coordinate searches,
queries of tables in a `specdb` database are managed
by a query *dict*.  The *dict* keys refer to the
table columns and the *dict* values are used to
perform the query.

The following table summarizes how different types
of values in a query *dict* are treated:

========== ============ =============================================
Value type Example      Interpretation
========== ============ =============================================
tuple      (1.,2.)      Interpreted as the min/max range for a query
int        1            Interpreted as a single value to match against
float      1.0          Interpreted as a single value to match against
str        '2'          Interpreted as a single value to match against
list       [2,32,256]   Interpreted as a list of values to match against
 ..         ..          Data type should match that of the table column
========== ============ =============================================

For *list* entries, the match is successful if the Table values
matches any of the items in the list.

Examples for Catalog Queries
----------------------------

While most queries of the database catalog will be
based on coordinates, one can further restrict the
query with a query *dict*.

Here are some simple examples for querying the :doc:`catalog`::

    qdict = {'zem': (1.,2.), 'STYPE': 'QSO'}

This will restrict to sources with redshift 1<=z<=2 and
source type 'QSO'.  Here is another::

    qdict = {'zem': (1.,2.), 'flag_group': 8}

Same cut on redshift and now restricting to sources
with the bitwise flag_group=8.

.. _bitwise-flags:

Bitwise Flags
=============

Extra handling is required to query bitwise flags
in a bitwise fashion.
The current approach is to append '-BITWISE' to the
keyword to indicate that the input value(s) is/are
to be matched in a bitwise manner.  In addition,
one appends '-AND' to require all of the bits are on
or '-OR' to require only one.  Here is an example::

    qdict = {'zem': (1.,2.), 'flag_group-BITWISE-OR': [2,4,32]}

This will match against the flag_group column with
bitwise logic.  That is, Table rows with any of the
input values turned on will match.

Note that the input values are the full values,
i.e.  2**5=32 or 2**9=512, not 5 or 9.
