#! /usr/bin/env python
# encoding: utf-8
#
# Copyright (C) 2011 Serge Monkewitz IPAC/Caltech
#

from __future__ import with_statement
import os
import sys
import traceback

from waflib import Build, Logs, Utils

def options(ctx):
    ctx.load('compiler_cxx CCfits boost armadillo')
    ctx.add_option('--debug', help='Include debug symbols and turn ' +
                                   'compiler optimizations off',
                   action='store_true', default=False, dest='debug')

def configure(ctx):
    ctx.load('compiler_cxx CCfits boost armadillo')
    ctx.env.append_value('CXXFLAGS', '-Wall')
    ctx.env.append_value('CXXFLAGS', '-Wextra')
    ctx.env.append_value('CXXFLAGS', '-std=c++11')
    ctx.env.append_value('CXXFLAGS', '-D__STDC_CONSTANT_MACROS')

    if ctx.options.debug:
        ctx.env.append_value('CXXFLAGS', '-g')
    else:
        ctx.env.append_value('CXXFLAGS', '-Ofast')
        ctx.env.append_value('CXXFLAGS', '-mtune=native')
        ctx.env.append_value('CXXFLAGS', '-march=native')
        ctx.env.append_value('CXXFLAGS', '-DNDEBUG')

    ctx.env.append_value('CXXFLAGS', '-std=c++11')

def build(ctx):
    cxx_sources=[
        'src/read_DRF.cxx',
        'src/Footprint/Footprint.cxx',
        'src/Footprint/compute_bounds.cxx',
        'src/Footprint/compute_correction.cxx',
        'src/Footprint/count_good_samples.cxx',
        'src/Footprint/get_response.cxx',
        'src/Footprint/generate_response.cxx',
        'src/Footprint/set_signals_to_sim_values.cxx',
        'src/Hires/compute_minimap.cxx',
        'src/Hires/iterate.cxx',
        'src/Hires/init.cxx',
        'src/Hires/write_output.cxx',
        'src/Hires/write_file.cxx',
        'src/Hires/spike_image.cxx',
        'src/Hires/start_image.cxx',
        'src/Hires/write_fits.cxx']
    headers=[
        'src/version.hxx',
        'src/Detector.hxx',
        'src/Gnomonic.hxx',
        'src/Hires.hxx',
        'src/Exception.hxx',
        'src/Sample.hxx',
        'src/Footprint.hxx']

    ctx.stlib(
         source=cxx_sources,
         target='hires',
         name='hires_st',
         install_path=os.path.join(ctx.env.PREFIX, 'lib'),
         use=['CCfits','boost','armadillo']
    )

    ctx.shlib(
        source=cxx_sources,
        target='hires',
        name='hires_sh',
        install_path=os.path.join(ctx.env.PREFIX, 'lib'),
        use=['CCfits','boost','armadillo']
    )

    ctx.install_files(
        os.path.join(ctx.env.PREFIX, 'include'),
        headers
        )
