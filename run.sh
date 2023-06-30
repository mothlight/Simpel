cd src
/usr/lib/jvm/default-java/bin/javac -d ../bin thud/simpel/SimpelModel.java
cd ../bin
directory=../example/
#readETFromFile=readET
readETFromFile=dontReadET
/usr/lib/jvm/default-java/bin/java thud.simpel.SimpelModel $directory $readETFromFile
