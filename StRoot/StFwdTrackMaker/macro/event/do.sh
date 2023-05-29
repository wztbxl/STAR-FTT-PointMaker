#!/bin/bash
input="event.lis"
iEvent=0
while IFS= read -r line
do
  echo "$line"
  echo root4star -b -q -l 'event/event_track.C( 500, "'"$line"'" )'
  root4star -b -q -l 'event/event_track.C( 500, "'"$line"'" )'
  iEvent=$((iEvent + 1))
  echo $iEvent
  mv FwdAna.root FwdAna_${iEvent}.root
done < "$input"