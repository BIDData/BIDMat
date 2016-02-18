set -e

#NDK_TOOLCHAIN_VERSION="snapdragonclang3.6"
#APP_ABI="armeabi-v7a"

cd jni/src/

ndk-build clean
ndk-build

cd ../..

#cp libs/$APP_ABI/*.so lib/
#rm -rf libs obj

#rm -rfv ../BIDMach_Android/app/src/main/libs/*
cp -rv ./jni/src/libs/* ../BIDMach_Android/app/src/main/libs/
cp -v ./BIDMat.jar ../BIDMach_Android/app/libs/
