#!/bin/bash

echo signal singlemuontop singleelectrontop singlemuonw singleelectrontop dimuon dielectron photon | xargs -n 1 -P 20 python makeLimitForest.py --region
