#region "copyright"

/*
    Copyright � 2016 - 2024 Stefan Berg <isbeorn86+NINA@googlemail.com> and the N.I.N.A. contributors

    This file is part of N.I.N.A. - Nighttime Imaging 'N' Astronomy.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#endregion "copyright"

using NINA.Core.Locale;
using NINA.Core.Utility;
using System;
using System.Globalization;
using System.Resources;
using System.Windows.Data;

namespace CanardConfit.NINA.BahtiFocus.Locale {

    public class BahtiLoc : BaseINPC, ILoc {
        private ResourceManager _locale;
        private CultureInfo _activeCulture;

        private static readonly Lazy<BahtiLoc> lazy = new(() => new BahtiLoc());

        private BahtiLoc() {
            _locale = new ResourceManager("CanardConfit.NINA.BahtiFocus.Locale.Locale", typeof(BahtiLoc).Assembly);
        }

        public void ReloadLocale(string culture) {
            using (MyStopWatch.Measure()) {
                try {
                    _activeCulture = new CultureInfo(culture);
                } catch (Exception ex) {
                    Logger.Error(ex);
                }
                RaiseAllPropertiesChanged();
            }
        }

        public static BahtiLoc Instance => lazy.Value;

        public string this[string key] {
            get {
                if (key == null) {
                    return string.Empty;
                }
                return _locale?.GetString(key, _activeCulture) ?? $"MISSING LABEL {key}";
            }
        }
    }

    public class BahtiLocExtension : Binding {

        public BahtiLocExtension(string name) : base($"[{name}]") {
            Mode = BindingMode.OneWay;
            Source = BahtiLoc.Instance;
        }
    }
}