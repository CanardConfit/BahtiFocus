using System.Reflection;
using System.Runtime.InteropServices;

// [MANDATORY] The following GUID is used as a unique identifier of the plugin. Generate a fresh one for your plugin!
[assembly: Guid("541625c3-50f8-4a60-86e1-57928fc0220a")]

// [MANDATORY] The assembly versioning
//Should be incremented for each new release build of a plugin
[assembly: AssemblyVersion("1.0.0.2")]
[assembly: AssemblyFileVersion("1.0.0.2")]

// [MANDATORY] The name of your plugin
[assembly: AssemblyTitle("BahtiFocus")]
// [MANDATORY] A short description of your plugin
[assembly: AssemblyDescription("BahtiFocus is a precision Bahtinov mask analyzer integrated into N.I.N.A., designed to help astronomers achieve perfect focus for their telescopes effortlessly.")]

// The following attributes are not required for the plugin per se, but are required by the official manifest metadata

// Your name
[assembly: AssemblyCompany("CanardConfit")]
// The product name that this plugin is part of
[assembly: AssemblyProduct("BahtiFocus")]
[assembly: AssemblyCopyright("Copyright © 2024 CanardConfit")]

// The minimum Version of N.I.N.A. that this plugin is compatible with
[assembly: AssemblyMetadata("MinimumApplicationVersion", "3.0.0.2017")]

// The license your plugin code is using
[assembly: AssemblyMetadata("License", "MPL-2.0")]
// The url to the license
[assembly: AssemblyMetadata("LicenseURL", "https://www.mozilla.org/en-US/MPL/2.0/")]
// The repository where your pluggin is hosted
[assembly: AssemblyMetadata("Repository", "https://github.com/CanardConfit/BahtiFocus")]

// The following attributes are optional for the official manifest metadata

//[Optional] Your plugin homepage URL - omit if not applicable
[assembly: AssemblyMetadata("Homepage", "https://github.com/CanardConfit/BahtiFocus")]

//[Optional] Common tags that quickly describe your plugin
[assembly: AssemblyMetadata("Tags", "Focus,Focus Analyser,Bahtinov Mask")]

//[Optional] A link that will show a log of all changes in between your plugin's versions
[assembly: AssemblyMetadata("ChangelogURL", "https://github.com/CanardConfit/BahtiFocus/CHANGELOG.md")]

//[Optional] The url to a featured logo that will be displayed in the plugin list next to the name
[assembly: AssemblyMetadata("FeaturedImageURL", "https://github.com/CanardConfit/BahtiFocus/blob/main/Images/logo.png?raw=true")]
//[Optional] A url to an example screenshot of your plugin in action
[assembly: AssemblyMetadata("ScreenshotURL", "https://github.com/CanardConfit/BahtiFocus/blob/main/Images/BahtiFocus_Working.png?raw=true")]
//[Optional] An additional url to an example screenshot of your plugin in action
[assembly: AssemblyMetadata("AltScreenshotURL", "https://github.com/CanardConfit/BahtiFocus/blob/main/Images/BahtiFocus_Start.png?raw=true")]
//[Optional] An in-depth description of your plugin
[assembly: AssemblyMetadata("LongDescription", @"*Prerequisites*

- [N.I.N.A 3.0](https://nighttime-imaging.eu/2024/03/18/n-i-n-a-3-0/)
- A connected camera (see [N.I.N.A documentation](https://nighttime-imaging.eu/docs/master/site/quickstart/equipment/))

*Usage*

- Install the plugin automatically with the N.I.N.A plugin manager or manually place the DLL file into `%localappdata%\NINA\Plugins`.
    - A new panel will appear in the `Imaging` tab called `Bahtinov Analyser`.
- Create or purchase a Bahtinov mask for your instrument (see [Bahtinov section](#bahtinov-mask)).
- Connect your camera to N.I.N.A.
- Choose, slew, and center a star in your camera's view.
    - If the star isn't bright enough, the lines created by the Bahtinov mask will not be visible.
- Place the Bahtinov mask on your instrument.
- Set the options in the `options` toolbar of the `Bahtinov Analyser` for exposure time, filters, etc.
    - **Note**: Information is taken from N.I.N.A options, such as focal length from `Options > Equipment`.
- Press the start button to begin the exposure.
- After the first photo is taken and appears on the right, you will see a rectangle on the right side of the screen.
    - The rectangle indicates the area analyzed on the left. You can move or zoom it using the controls at the bottom to center it on the star.
- The analysis will continue automatically, and you can adjust the focus of your instrument to find the critical focus.
- Once the critical focus is achieved, you can stop the process and remove the Bahtinov mask!
")]

// Setting ComVisible to false makes the types in this assembly not visible
// to COM components.  If you need to access a type in this assembly from
// COM, set the ComVisible attribute to true on that type.
[assembly: ComVisible(false)]

// [Unused]
[assembly: AssemblyConfiguration("")]
// [Unused]
[assembly: AssemblyTrademark("")]
// [Unused]
[assembly: AssemblyCulture("")]