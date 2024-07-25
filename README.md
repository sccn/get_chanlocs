![Screenshot 2024-07-25 at 14 57 29](https://github.com/user-attachments/assets/fc3850a2-09d2-4341-aecf-a41fc39e0758)

# The get_chanlocs EEGLAB plugin
The get_chanlocs EEGLAB plug-in is built on functions in FieldTrip to locate 3-D electrode positions from a 3-D scanned head image. Robert Oostenveld, originator of the FieldTrip toolbox, alerted us in 2017 that he and his students in Nijmegen had put functions into FieldTrip to compute positions of scalp electrodes from the recorded 3-D images for one 3-D camera, the Structure scanner mounted to an Apple iPad. (Read Homölle and Oostenveld (2019) and notes on the incorporated FieldTrip functions). We at SCCN have created an EEGLAB plug-in extension, get_chanlocs, to ease the process of digitizing the positions of the electrodes from the acquired 3-D and entering them into the EEG.chanlocs data structure for use with other EEGLAB (plotting and source localization) functions that require electrode position information.

See the [wiki](https://github.com/sccn/get_chanlocs/wiki) for more information or the other menu for this plugin if you are on the EEGLAB web page.
