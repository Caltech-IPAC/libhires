#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    # Find libxml2
    if conf.options.libxml2_dir:
        if not conf.options.libxml2_incdir:
            conf.options.libxml2_incdir=conf.options.libxml2_dir + "/include"
        if not conf.options.libxml2_libdir:
            conf.options.libxml2_libdir=conf.options.libxml2_dir + "/lib"

    libxml2_incdir=[None,'/usr/include/libxml2']
    if conf.options.libxml2_incdir:
        libxml2_incdir=[conf.options.libxml2_incdir]
        
    libxml2_libdir=[]
    if conf.options.libxml2_libdir:
        libxml2_libdir=[conf.options.libxml2_libdir]

    libxml2_libs=['xml2']
    if conf.options.libxml2_libs:
        libxml2_libs=conf.options.libxml2_libs.split()

    found_libxml2=False
    for incdir in libxml2_incdir:
        try:
            conf.check_cc(msg="Checking for libxml2 using include = " + str(incdir),
                          header_name='libxml/parser.h',
                          includes=[incdir],
                          uselib_store='libxml2',
                          libpath=libxml2_libdir,
                          rpath=libxml2_libdir,
                          lib=libxml2_libs)
        except conf.errors.ConfigurationError:
            continue
        else:
            found_libxml2=True
            break
    if not found_libxml2:
        conf.fatal("Could not find libxml2 libraries")
        
def options(opt):
    libxml2=opt.add_option_group('Libxml2 Options')
    libxml2.add_option('--libxml2-dir',
                   help='Base directory where libxml2 is installed')
    libxml2.add_option('--libxml2-incdir',
                   help='Directory where libxml2 include files are installed')
    libxml2.add_option('--libxml2-libdir',
                   help='Directory where libxml2 library files are installed')
    libxml2.add_option('--libxml2-libs',
                   help='Names of the libxml2 libraries without prefix or suffix\n'
                   '(e.g. "xml2")')
