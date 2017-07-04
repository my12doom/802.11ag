# 802.11ag
yet another software defined radio 802.11a/g tranciver

----------
This is a project for Coded-OFDM radio level STUDY, from scratch.

During initial research, I found a lot of useless bullshit articals on search engine.  
To avoid being theoretical talking or pure simulation, the first milestone is to receive OFDM signals sent by real life commercial wifi cards.

The final goal is to understand the philosophy behind wifi standards design, and then design a new baseband that improves its performance and/or prices for new application circustance.
(Yes, I started this because i'm angry about shitty wifi "HD Streaming" on most drones.) 


----------
# Milestones
of course all things down below are to be done by SDR, not wifi cards.

+ Receive correct data frames from commercial wifi cards. (It works!)
+ Receive live video stream from "wifibroadcast" (optimization)
 (see [https://befinitiv.wordpress.com/wifibroadcast-analog-like-transmission-of-live-video-data/](https://befinitiv.wordpress.com/wifibroadcast-analog-like-transmission-of-live-video-data/ "https://befinitiv.wordpress.com/wifibroadcast-analog-like-transmission-of-live-video-data/") for wifibroadcast)
+ Receive wifibroadcast on Android / embedded device (more optimization and cross platform). 
+ Understand tradeoffs in wifi standards, improve them and implement a better data link (goal).
+ FPGA/Hardware implementations (dream)

----------
# Coding
+ Standard C++ or assembly, no external libraries.