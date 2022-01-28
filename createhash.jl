# Create a hash of the log for an untainted webapp container
# It is important that the log is not stateful, i.e. that the 
# log output is always the same after the container is fully initialized. 

using SHA

file = "swarm_virtualdma"    # hashfile
ID = "3733099a9c06"          # ID of an untainted fully executed app
hashlog(ID) = (read(`docker logs $(ID)`) |> sha256 |> bytes2hex)

open(f -> write(f, hashlog(ID)), file, "w")
