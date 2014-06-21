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
        'src/cxx/read_all_DRF_files/read_all_DRF_files.cxx',
        'src/cxx/read_all_DRF_files/read_all_DRF_planck.cxx',
        'src/cxx/read_all_IN_files/read_all_IN_files.cxx',
        'src/cxx/read_all_IN_files/read_one_IN_planck.cxx',
        'src/cxx/Footprint/Footprint.cxx',
        'src/cxx/Footprint/compute_bounds.cxx',
        'src/cxx/Footprint/compute_correction.cxx',
        'src/cxx/Footprint/count_good_samples.cxx',
        'src/cxx/Footprint/get_response.cxx',
        'src/cxx/Footprint/generate_response.cxx',
        'src/cxx/Footprint/set_fluxes_to_sim_values.cxx',
        'src/cxx/Footprint/compute_minimap.cxx',
        'src/cxx/Hires_Parameters/set_default_values.cxx',
        'src/cxx/Hires_Parameters/parse_command_line.cxx',
        'src/cxx/Hires_Parameters/dump_params.cxx',
        'src/cxx/Hires/iterate.cxx',
        'src/cxx/Hires/write_output.cxx',
        'src/cxx/Hires/ostream_operator.cxx',
        'src/cxx/Hires/spike_image.cxx',
        'src/cxx/Hires/start_image.cxx',
        'src/cxx/Hires/write_fits.cxx']
    headers=[
        'src/cxx/version.hxx',
        'src/cxx/Detector.hxx',
        'src/cxx/Gnomonic.hxx',
        'src/cxx/Hires_Parameters.hxx',
        'src/cxx/Hires.hxx',
        'src/cxx/Exception.hxx',
        'src/cxx/Sample.hxx',
        'src/cxx/Footprint.hxx']

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

#    ctx.program(
#        source=['src/cxx/main.cxx'],
#        target='hires',
#        name='hires_bin',
#        lib=['cfitsio'],
#        install_path=os.path.join(ctx.env.PREFIX, 'bin'),
#        use=['CCfits','boost','hires_sh','armadillo']
#    )

    ctx.install_files(
        os.path.join(ctx.env.PREFIX, 'include'),
        headers
        )
