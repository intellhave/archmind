#Detect the number of cpus
import os
num_cpu = int(os.environ.get('NUM_CPU',2))
SetOption('num_jobs', num_cpu+1)

#Build options
vars = Variables()
vars.Add('puthonlib', 'Set to 1 to build a python lib', 0)
vars.Add('debug', 'Set to 1 to build with multi precision', 0)

env = Environment( 
          variables = vars,
          SHLIBPREFIX = '',	#do not prefix the library with lib
	  LIBS = Split('boost_python python2.6'),
	  LIBPATH = Split('#/lib /usr/lib'),
          CPPPATH = Split('#/include /usr/include /usr/include/python2.6'),
        )

#Generate help for the options
Help(vars.GenerateHelpText(env))

pythonlib = ARGUMENTS.get('pythonlib',0)    
debug = ARGUMENTS.get('debug',0)

Export('pythonlib')

conf = Configure(env)
env = conf.Finish()

if debug:
    env.Append(CPPFLAGS = Split('-Wall -g '))
else:
    env.Append(CPPFLAGS = Split('-Wall -DNDEBUG -O3 -msse -msse2'))

#add modules from respective files   
objs = []

objs.append(SConscript('src/Python/SConscript', exports='env'))

if pythonlib:
    env.SharedLibrary(target = './Lib/melina',source = objs)
else:
    objs.append(SConscript('src/App/SConscript', exports='env'))
    env.Program(target = './Release/app',source = objs)


