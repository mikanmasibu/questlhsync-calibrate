# QuestLHSync Calibrate

PC SteamVR driver that keeps lighthouse trackers lined up with a Quest or a Steam Frame, and adds a one-time calibration for the tilt and offset that the camera fit leaves behind.

It is built from two projects:

- **[QuestLHSync](https://github.com/CreoleVR/QuestLHSync)** by CreoleVR (MIT). This is that driver. The headset's tracking cameras watch the base stations' laser flashes and solve the alignment while you play. Two additions are on top: a settled alignment ignores a sudden ~180 degree swap of two base stations, and a one-time calibration can take out a small remaining tilt and shift.
- **[OpenVR-SpaceSync](https://github.com/shinyflvre/OpenVR-SpaceSync)** by shinyflvre (AGPL-3.0). SpaceSync, itself based on [OpenVR-SpaceOverride](https://github.com/Nyabsi/OpenVR-SpaceOverride) and [OpenVR-SpaceCalibrator](https://github.com/pushrax/OpenVR-SpaceCalibrator), calibrates by holding a lighthouse tracker on the headset and looking left, right, up and down. This build uses that same capture. The solver here was written for QuestLHSync. It is not a copy of SpaceSync's source.

Do not run SpaceSync, Space Calibrator, or another tool that also moves lighthouse devices. Two of them fight.

The code in this repository is MIT, the same license as QuestLHSync. See [LICENSE](LICENSE).

## Install

**PC.** Download `QuestLHSync-Calibrate-Installer.exe` from [Releases](https://github.com/mikanmasibu/questlhsync-calibrate/releases), run it, then start SteamVR. It installs only the SteamVR driver, to `%LOCALAPPDATA%\QuestLHSync\questlhsync`, and registers it. SteamVR can be open; the installer closes it.

**Quest or Steam Frame.** The headset software is not in this installer. Install it from the original QuestLHSync release, then use this PC driver with it:

[github.com/CreoleVR/QuestLHSync/releases](https://github.com/CreoleVR/QuestLHSync/releases)

- Quest (Pro, 3 or 3S, rooted with Magisk): `QuestLHSync-quest-module-<version>.zip`
- Steam Frame: `QuestLHSync-frame-installer.flatpak` or `QuestLHSync-frame-module-<version>.tar.gz`

Follow the original install notes for the headset. The PC side of that release is not needed; this installer replaces it.

## One-time calibration

Wait until the dashboard says the alignment is locked and settled. Hold a tracker or a controller firmly on the headset, click **Calibrate**, and look left, right, up and down over about 12 seconds. The correction (degrees of tilt, and centimetres at the head) is saved and applied on top of the live camera alignment. **Clear** removes it.

A correction larger than 8 degrees or 8 cm is rejected, as is a device that was not held on the headset.

The PC driver follows QuestLHSync 1.16. Use the 1.16 headset package from the original release. A wired headset (Index, Vive, and other lighthouse headsets) turns this driver off. Base stations are held at the average of SteamVR's measurements (`steadyStations`, on by default). The installer sets SteamVR's `activateMultipleDrivers` to true so this driver loads beside the lighthouse driver.

## Building

Visual Studio 2022 with C++:

```
build.bat
```

The installer is `out\QuestLHSync-Calibrate-Installer.exe`. It contains the driver just built.
