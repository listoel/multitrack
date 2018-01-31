********
Tutorial
********

The basic building blocks of a ``multitrack`` simulation are a ring and
an extraction point, defined through the :class:`~multitrack.inputs.Ring`
and :class:`~multitrack.inputs.Extraction` classes, together with a
particle distribution defined with :func:`~multitrack.inputs.get_init()`.
These come together in the :func:`~multitrack.track.track()` function.
