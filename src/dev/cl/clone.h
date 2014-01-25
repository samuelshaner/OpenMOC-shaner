#include "../DeviceMaterial.h"
#include "../DeviceTrack.h"

void clCloneMaterial(Material* material_h, dev_material* material_d);
void clCloneTrackOnGPU(Track* track_h, dev_track* track_d);
