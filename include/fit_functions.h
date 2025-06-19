#pragma once

TF1 *bkg = new TF1("bkg", "[0]*TMath::Power(x,[1])*(1-TMath::Erf((x-[2])/[3]))", 0, 150);
TF1 *sig = new TF1("sig", "[0]*TMath::Power(x,[1])*(1-TMath::Erf((x-[2])/[3]))+gaus(4)", 0, 150);
bkg->SetParameters(1.86764, -1.11647, 100.758, 11.8833);

/*

p0                        =      1.86764   +/-   1.67246e-07
p1                        =     -1.11647   +/-   2.62015e-08
p2                        =      100.758   +/-   0.492804
p3                        =      11.8833   +/-   0.807918

*/

double _h_find_center(double *variable_s, double *parameter_s)
{
    //  --- Variables
    auto _x = variable_s[0];
    auto _y = variable_s[1];

    //  --- Signal Parameters
    auto _center_x = parameter_s[0];
    auto _center_y = parameter_s[1];
    auto _radius = parameter_s[2];
    auto _sigma = parameter_s[3];
    auto _amplitude = parameter_s[4];

    //  --- Background Parameters
    auto _bkg_amplitude = parameter_s[5];
    auto _bkg_power = parameter_s[6];
    auto _bkg_drop = parameter_s[7];
    auto _bkg_drop_sigma = parameter_s[8];

    //  --- Radial variables
    auto _R = TMath::Sqrt((_x - _center_x) * (_x - _center_x) + (_y - _center_y) * (_y - _center_y));
    auto _R_no_center = TMath::Sqrt((_x) * (_x) + (_y) * (_y));
    auto _phi = TMath::ATan2((_y - _center_y), (_x - _center_x));

    if (_R_no_center > 100 || _R_no_center < 30)
        return 0;

    //  --- Function
    auto _signal = _amplitude * TMath::Gaus(_R, _radius, _sigma, kTRUE);
    auto _background = _bkg_amplitude * TMath::Power(_R, _bkg_power) * (1 - TMath::Erf((_R_no_center - _bkg_drop) / _bkg_drop_sigma));
    return _signal + _background;
}
double _h_find_center_2(double *variable_s, double *parameter_s)
{

    //  --- Variables
    auto _x = variable_s[0];
    auto _y = variable_s[1];

    //  --- Signal Parameters
    auto _center_x = parameter_s[0];
    auto _center_y = parameter_s[1];
    auto _radius = parameter_s[9];
    auto _sigma = parameter_s[3];
    auto _amplitude = parameter_s[10];

    //  --- Radial variables
    auto _R = TMath::Sqrt((_x - _center_x) * (_x - _center_x) + (_y - _center_y) * (_y - _center_y));

    //  --- Function
    auto _signal = _amplitude * TMath::Gaus(_R, _radius, _sigma, kTRUE);
    auto _background = _h_find_center(variable_s, parameter_s);
    return _signal + _background;
}
double _h_find_center_3(double *variable_s, double *parameter_s)
{

    //  --- Variables
    auto _x = variable_s[0];
    auto _y = variable_s[1];

    //  --- Signal Parameters
    auto _center_x = parameter_s[0];
    auto _center_y = parameter_s[1];
    auto _radius = parameter_s[11];
    auto _sigma = parameter_s[3];
    auto _amplitude = parameter_s[12];

    //  --- Radial variables
    auto _R = TMath::Sqrt((_x - _center_x) * (_x - _center_x) + (_y - _center_y) * (_y - _center_y));

    //  --- Function
    auto _signal = _amplitude * TMath::Gaus(_R, _radius, _sigma, kTRUE);
    auto _background = _h_find_center_2(variable_s, parameter_s);
    return _signal + _background;
}

//  TODO
std::array<float, 2>
eval_with_errors(TGraphErrors *gTarget, float _xtarget)
{
    auto iPnt = 0;
    auto current_x = 1.;
    auto previous_x = 1.;
    float interpolated_y = -1.;
    float interpolated_ey = -1.;
    for (iPnt = 1; iPnt < gTarget->GetN() - 1; iPnt++)
    {
        current_x = gTarget->GetPointX(iPnt);
        previous_x = gTarget->GetPointX(iPnt - 1);
        if ((current_x >= _xtarget) && (previous_x <= _xtarget))
        {
            break;
        }
    }
    if (iPnt == gTarget->GetN() - 1)
    {
        return {interpolated_y, interpolated_ey};
    }
    auto coeff_y1 = (_xtarget - previous_x) / (current_x - previous_x);
    auto coeff_y0 = (current_x - _xtarget) / (current_x - previous_x);
    interpolated_y = coeff_y1 * gTarget->GetPointY(iPnt) + coeff_y0 * gTarget->GetPointY(iPnt - 1);
    interpolated_ey = sqrt((coeff_y0 * gTarget->GetErrorY(iPnt)) * (coeff_y0 * gTarget->GetErrorY(iPnt)) + (coeff_y1 * gTarget->GetErrorY(iPnt - 1)) * (coeff_y1 * gTarget->GetErrorY(iPnt - 1)));
    return {interpolated_y, interpolated_ey};
}
std::array<float, 2>
eval_with_errors(TGraph *gTarget, float _xtarget)
{
    auto iPnt = 0;
    auto current_x = 1.;
    auto previous_x = 1.;
    float interpolated_y = -1.;
    float interpolated_ey = -1.;
    for (iPnt = 1; iPnt < gTarget->GetN() - 1; iPnt++)
    {
        current_x = gTarget->GetPointX(iPnt);
        previous_x = gTarget->GetPointX(iPnt - 1);
        if ((current_x >= _xtarget) && (previous_x <= _xtarget))
        {
            break;
        }
    }
    if (iPnt == gTarget->GetN() - 1)
    {
        return {interpolated_y, interpolated_ey};
    }
    auto coeff_y1 = (_xtarget - previous_x) / (current_x - previous_x);
    auto coeff_y0 = (current_x - _xtarget) / (current_x - previous_x);
    interpolated_y = coeff_y1 * gTarget->GetPointY(iPnt) + coeff_y0 * gTarget->GetPointY(iPnt - 1);
    interpolated_ey = 0.;
    return {interpolated_y, interpolated_ey};
}

int which_zone(std::pair<std::array<float, 2>, std::vector<std::array<float, 2>>> zones_delimiter, float target_x, float target_y)
{
    //  Calculate radius
    auto target_radius = sqrt((target_x - zones_delimiter.first[0]) * (target_x - zones_delimiter.first[0]) + (target_y - zones_delimiter.first[1]) * (target_y - zones_delimiter.first[1]));
    auto zone = 0;
    for (auto i = 0; i < zones_delimiter.second.size(); i++)
    {
        if (target_radius >= zones_delimiter.second[i][0] && target_radius < zones_delimiter.second[i][1])
        {
            zone = i + 1;
            break;
        }
    }
    return zone;
}

float which_sensor(std::array<float, 2> position)
{
    if (position[1] > 28)
        return 0;
    if (position[1] < -28)
        return 1;
    if (position[0] < -28)
        return 1;
    return 0;
}