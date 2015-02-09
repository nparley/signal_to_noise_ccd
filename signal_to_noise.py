"""
Flask signal to noise calculator. Allows you to calculate the signal-to-noise ratio for a star
observed with a telescope and CCD. You can vary the parameters and find out what sort of results
to expect. Code adapted from original perl code
here: http://stupendous.rit.edu/richmond/signal.shtml
"""
import flask
import flask.ext.wtf
from wtforms.fields import FloatField, SelectField, SubmitField
from wtforms.validators import DataRequired, NumberRange
import functions

app = flask.Flask(__name__)
app.debug = True
app.config['SECRET_KEY'] = 'put secret key here'


class InputForm(flask.ext.wtf.Form):
    filter_name = SelectField('Filter. This calculator uses Bessell filters, which are '
                              'close to Johnson UBV, Cousins RI', choices=[('none', 'none'),
                                                                           ('U', 'U'), ('B', 'B'),
                                                                           ('V', 'V'), ('R', 'R'),
                                                                           ('I', 'I')])
    mag_start = FloatField('Start Magnitude',  validators=[DataRequired()], default=10.0)
    mag_end = FloatField('End Magnitude', validators=[DataRequired()], default=11.0)
    dmag = FloatField('Magnitude steps', validators=[DataRequired()], default=0.5)
    tel_diam = FloatField('Telescope diameter (cm)', validators=[DataRequired()], default=10.0)
    qe = FloatField('Overall QE (0.0-1.0), e.g. 0.3 for cheap CCD / 0.7 expensive CCD',
                    validators=[DataRequired(), NumberRange(min=0.0, max=1.0)], default=0.5)
    pixsize = FloatField('Pixel size (arcsec/pixel)',
                         validators=[DataRequired()], default=3.5)
    readnoise = FloatField('CCD readout noise (ADU)',
                           validators=[DataRequired()], default=5.0)
    skymag = FloatField('Sky brightness (mag per sq. arcsec) e.g. Mountains 21, suburbs 17',
                        validators=[DataRequired()], default=21)
    airmass = FloatField('Airmass', validators=[DataRequired()], default=1.2)
    exptime = FloatField('Exposure time (seconds)', validators=[DataRequired()], default=30.0)
    fwhm = FloatField('FWHM (arcsec)', validators=[DataRequired()], default=8.0)
    aper_rad = FloatField('Radius for photometry area over which signal is measured (arcsec)',
                          validators=[DataRequired()], default=10.0)
    obs_diam = FloatField('Obscuration diameter (cm)', default=0.0)
    gain = FloatField('CCD gain electrons / ADU', validators=[DataRequired()], default=1.0)
    sat_level = FloatField('Saturation level (ADU)', validators=[DataRequired()], default=65000)
    submit = SubmitField('Calculate Results')


@app.route('/', methods=('GET', 'POST'))
def create_input_form():
    form = InputForm()
    if form.validate_on_submit():
        try:
            output = functions.signal_to_noise(filter_name=form.filter_name.data,
                                               mag_start=form.mag_start.data,
                                               mag_end=form.mag_end.data, dmag=form.dmag.data,
                                               tel_diam=form.tel_diam.data, qe=form.qe.data,
                                               readnoise=form.readnoise.data,
                                               pixsize=form.pixsize.data,
                                               skymag=form.skymag.data, airmass=form.airmass.data,
                                               gain=form.gain.data, obs_diam=form.obs_diam.data,
                                               exptime=form.exptime.data, fwhm=form.fwhm.data,
                                               sat_level=form.sat_level.data,
                                               aper_rad=form.aper_rad.data)

            return flask.render_template('results.html', results=output)

        except ValueError as e:
            return flask.render_template('error.html', error_message=e)
    else:
        return flask.render_template('input.html', form=form)


if __name__ == '__main__':
    app.run()
