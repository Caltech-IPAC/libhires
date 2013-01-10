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
    ctx.load('compiler_cxx')
    ctx.add_option('--debug', help='Include debug symbols and turn ' +
                                   'compiler optimizations off',
                   action='store_true', default=False, dest='debug')

    ccfits=ctx.add_option_group('CCfits Options')
    ccfits.add_option('--ccfits-dir',
                   help='Base directory where ccfits is installed')
    ccfits.add_option('--ccfits-incdir',
                   help='Directory where ccfits include files are installed')
    ccfits.add_option('--ccfits-libdir',
                   help='Directory where ccfits library files are installed')
    ccfits.add_option('--ccfits-libs',
                   help='Names of the ccfits libraries without prefix or suffix\n'
                   '(e.g. "CCFITS"')

    log4cxx=ctx.add_option_group('log4cxx Options')
    log4cxx.add_option('--log4cxx-dir',
                   help='Base directory where log4cxx is installed')
    log4cxx.add_option('--log4cxx-incdir',
                   help='Directory where log4cxx include files are installed')
    log4cxx.add_option('--log4cxx-libdir',
                   help='Directory where log4cxx library files are installed')
    log4cxx.add_option('--log4cxx-libs',
                   help='Names of the log4cxx libraries without prefix or suffix\n'
                   '(e.g. "log4cxx"')

    boost=ctx.add_option_group('boost Options')
    boost.add_option('--boost-dir',
                   help='Base directory where boost is installed')
    boost.add_option('--boost-incdir',
                   help='Directory where boost include files are installed')
    boost.add_option('--boost-libdir',
                   help='Directory where boost library files are installed')
    boost.add_option('--boost-libs',
                   help='Names of the boost libraries without prefix or suffix\n'
                   '(e.g. "boost_filesystem boost_system"')

def configure(ctx):
    ctx.load('compiler_cxx')
    ctx.env.append_value('CXXFLAGS', '-Wall')
    ctx.env.append_value('CXXFLAGS', '-Wextra')
    ctx.env.append_value('CXXFLAGS', '-std=c++11')
    ctx.env.append_value('CXXFLAGS', '-D__STDC_CONSTANT_MACROS')

    # Find CCFITS
    if ctx.options.ccfits_dir:
        if not ctx.options.ccfits_incdir:
            ctx.options.ccfits_incdir=ctx.options.ccfits_dir + "/include"
        if not ctx.options.ccfits_libdir:
            ctx.options.ccfits_libdir=ctx.options.ccfits_dir + "/lib"
    frag="#include <CCfits>\n" + 'int main()\n' \
        + "{CCfits::FITS::setVerboseMode(true);}\n"
    if ctx.options.ccfits_incdir:
        ccfits_inc=ctx.options.ccfits_incdir
    else:
        ccfits_inc='/usr/include/CCfits'
    if ctx.options.ccfits_libs:
        ccfits_libs=ctx.options.ccfits_libs
    else:
        ccfits_libs="CCfits"

    ctx.check_cxx(msg="Checking for CCfits",
                  fragment=frag,
                  includes=[ccfits_inc], uselib_store='ccfits',
                  libpath=[ctx.options.ccfits_libdir],
                  rpath=[ctx.options.ccfits_libdir],
                  lib=[ccfits_libs])


    # Find LOG4CXX
    if ctx.options.log4cxx_dir:
        if not ctx.options.log4cxx_incdir:
            ctx.options.log4cxx_incdir=ctx.options.log4cxx_dir + "/include"
        if not ctx.options.log4cxx_libdir:
            ctx.options.log4cxx_libdir=ctx.options.log4cxx_dir + "/lib"
    frag="#include <log4cxx/logger.h>\n" + 'int main()\n' \
        + "{log4cxx::Logger::getRootLogger();}\n"
    if ctx.options.log4cxx_incdir:
        log4cxx_inc=ctx.options.log4cxx_incdir
    else:
        log4cxx_inc='/usr/include'
    if ctx.options.log4cxx_libs:
        log4cxx_libs=ctx.options.log4cxx_libs
    else:
        log4cxx_libs="log4cxx"

    ctx.check_cxx(msg="Checking for Log4cxx",
                  fragment=frag,
                  includes=[log4cxx_inc], uselib_store='log4cxx',
                  libpath=[ctx.options.log4cxx_libdir],
                  rpath=[ctx.options.log4cxx_libdir],
                  lib=[log4cxx_libs])

    # Find Boost
    if ctx.options.boost_dir:
        if not ctx.options.boost_incdir:
            ctx.options.boost_incdir=ctx.options.boost_dir + "/include"
        if not ctx.options.boost_libdir:
            ctx.options.boost_libdir=ctx.options.boost_dir + "/lib"
    frag="#include <boost/filesystem.hpp>\n" + 'int main()\n' \
        + "{boost::filesystem::path();}\n"
    if ctx.options.boost_incdir:
        boost_inc=ctx.options.boost_incdir
    else:
        boost_inc='/usr/include'
    if ctx.options.boost_libs:
        boost_libs=[ctx.options.boost_libs]
    else:
        boost_libs=["boost_filesystem","boost_system"]

    ctx.check_cxx(msg="Checking for Boost",
                  fragment=frag,
                  includes=[boost_inc], uselib_store='boost',
                  libpath=[ctx.options.boost_libdir],
                  rpath=[ctx.options.boost_libdir],
                  lib=boost_libs)

    if ctx.options.debug:
        ctx.env.append_value('CXXFLAGS', '-g')
    else:
        ctx.env.append_value('CXXFLAGS', '-O2')
        ctx.env.append_value('CXXFLAGS', '-mtune=native')
        ctx.env.append_value('CXXFLAGS', '-march=native')
        ctx.env.append_value('CXXFLAGS', '-DNDEBUG')

    ctx.env.append_value('CXXFLAGS', '-std=c++11')

    ctx.env.append_value('CPPFLAGS', '-I/usr/include/eigen3')

    ctx.env.append_value('LINKFLAGS', '-larmadillo')

def build(ctx):
    # tree index generator
    ctx.program(
        source=[
            'src/cxx/main.cxx',
            'src/cxx/compute_correction.cxx',
            'src/cxx/make_start_image.cxx',
            'src/cxx/read_all_DRF_files/read_all_DRF_files.cxx',
            'src/cxx/read_all_DRF_files/read_all_DRF_planck.cxx',
            'src/cxx/read_all_IN_files/read_all_IN_files.cxx',
            'src/cxx/read_all_IN_files/read_one_IN_planck.cxx',
            'src/cxx/Footprint/Footprint.cxx',
            'src/cxx/Footprint/compute_bounds.cxx',
            'src/cxx/Footprint/count_good_samples.cxx',
            'src/cxx/Footprint/get_response.cxx',
            'src/cxx/Footprint/generate_response.cxx',
            'src/cxx/Params/Params.cxx',
            'src/cxx/Params/ostream_operator.cxx',
            'src/cxx/write_fits.cxx'],
        target='hires',
        name='hires',
        install_path=os.path.join(ctx.env.PREFIX, 'bin'),
        use=['ccfits','log4cxx','boost']
    )

