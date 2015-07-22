#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    conf.load(['armadillo','libxml2'])
    # Find MLPACK
    if conf.options.mlpack_dir:
        if not conf.options.mlpack_incdir:
            conf.options.mlpack_incdir=conf.options.mlpack_dir + "/include"
        if not conf.options.mlpack_libdir:
            conf.options.mlpack_libdir=conf.options.mlpack_dir + "/lib"

    if conf.options.mlpack_incdir:
        mlpack_incdir=[conf.options.mlpack_incdir]
    else:
        mlpack_incdir=[]
    if conf.options.mlpack_libdir:
        mlpack_libdir=[conf.options.mlpack_libdir]
    else:
        mlpack_libdir=[]

    if conf.options.mlpack_libs:
        mlpack_libs=conf.options.mlpack_libs.split()
    else:
        mlpack_libs=['mlpack']

    conf.check_cxx(msg="Checking for mlpack",
                   header_name='mlpack/core.hpp',
                   includes=mlpack_incdir,
                   uselib_store='mlpack',
                   libpath=mlpack_libdir,
                   rpath=mlpack_libdir,
                   lib=mlpack_libs,
                   use=['armadillo','libxml2'])

def options(opt):
    mlpack=opt.add_option_group('Mlpack Options')
    mlpack.add_option('--mlpack-dir',
                   help='Base directory where mlpack is installed')
    mlpack.add_option('--mlpack-incdir',
                   help='Directory where mlpack include files are installed')
    mlpack.add_option('--mlpack-libdir',
                   help='Directory where mlpack library files are installed')
    mlpack.add_option('--mlpack-libs',
                   help='Names of the mlpack libraries without prefix or suffix\n'
                   '(e.g. "mlpack")')
    opt.load(['armadillo','libxml2'])
