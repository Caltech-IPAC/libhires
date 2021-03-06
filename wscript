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
    ctx.load(['compiler_cxx','cxx11','cfitsio','CCfits','boost',
              'mlpack','armadillo','libxml2'])
    ctx.add_option('--debug', help='Include debug symbols and turn ' +
                                   'compiler optimizations off',
                   action='store_true', default=False, dest='debug')

def configure(ctx):
    ctx.load(['compiler_cxx','cxx11','cfitsio','CCfits','boost',
              'mlpack','armadillo','libxml2'])
    ctx.check_boost('filesystem system')
    ctx.env.append_value('CXXFLAGS', '-Wall')
    ctx.env.append_value('CXXFLAGS', '-Wextra')
    if ctx.options.debug:
        ctx.env.append_value('CXXFLAGS', '-g')
    else:
        ctx.env.append_value('CXXFLAGS', '-Ofast')
        ctx.env.append_value('CXXFLAGS', '-DNDEBUG')
        ctx.env.append_value('CXXFLAGS', '-DARMA_NO_DEBUG')

def build(ctx):
    cxx_sources=[
        'src/Binned_Data.cxx',
        'src/Hires/compute_minimap.cxx',
        'src/Hires/compute_mcm.cxx',
        'src/Hires/compute_elastic_net.cxx',
        'src/Hires/compute_tikhonov.cxx',
        'src/Hires/compute_response_function.cxx',
        'src/Hires/Hires.cxx',
        'src/Hires/write_file.cxx',
        'src/Hires/write_fits.cxx']
    headers=[
        'src/version.hxx',
        'src/Gnomonic.hxx',
        'src/Hires.hxx',
        'src/Binned_Data.hxx',
        'src/Exception.hxx',
        'src/Sample.hxx']

    ctx.stlib(
         source=cxx_sources,
         target='hires',
         name='hires_st',
         install_path=os.path.join(ctx.env.PREFIX, 'lib'),
         use=['cxx11','cfitsio','CCfits','BOOST','mlpack','armadillo',
              'libxml2']
    )

    ctx.shlib(
        source=cxx_sources,
        target='hires',
        name='hires_sh',
        install_path=os.path.join(ctx.env.PREFIX, 'lib'),
        use=['cxx11','cfitsio','CCfits','BOOST','mlpack','armadillo',
             'libxml2']
    )

    ctx.program(source=['test/deconvolve.cxx'],
                target='deconvolve',
                installpath=None,
                use=['cxx11','hires_st','BOOST']
            )
                
    ctx.install_files(
        os.path.join(ctx.env.PREFIX, 'include'),
        headers
        )
