# otf2ttf
Turn OTF to TTF using [caryll/otfcc](https://github.com/caryll/otfcc "caryll/otfcc").

Migrate from [caryll/otfcc-cubic2quad](https://github.com/caryll/otfcc-cubic2quad "caryll/otfcc-cubic2quad") to java version, and fixed the problem that can not convert large OTF files.

## Installation
1. Install [JDK8](http://www.oracle.com/technetwork/java/javase/downloads/index.html) and [Apache Maven](http://maven.apache.org).
1. Download [caryll/otfcc](https://github.com/caryll/otfcc "caryll/otfcc") prebuilt binaries and unpacking.
1. Build otf2ttf `git clone https://github.com/icaoweiwei/otf2ttf.git && cd otf2ttf && mvn clean install`.


## Usage
1. Copy `target/otf2ttf.jar` and `build.bat` to `otfcc prebuilt binaries` folder.
1. Run `build.bat <folder path>` to convert all OTF files in folder path.
