__author__ = 'Alessandro Orchini'

import numpy as np
import time


def sinefit(*args):
    """Fit a sinusoidal function"""
    # params = sinefit(yin,t,f)
    # [params,yest] = sinefit(yin,t,f)
    # [params,yest] = sinefit(yin,t,f,verbose)
    # [params,yest] = sinefit(yin,t,f,verbose,plotflag)
    # [params,yest] = sinefit(yin,t,f,verbose,plotflag,searchflag)
    # [params,yest,yres] = sinefit(yin,t,...)
    # [params,yest,yres,rmserr] = sinefit(yin,t,...)
    #
    # Least squares sinusoid fit - optimization toolbox not needed
    #
    # IEEE Standard for Digitizing Waveform Recorders (IEEE Std 1057):
    # Algorithm for three-parameter (known frequency) and four-parameter
    # (general use) least squares fit to sinewave data using matrix operations.
    #
    # INPUTS:   -yin:        Input signal to be fitted (a vector)
    # -t:          Time vector (same length as yin)
    #           -f:          Signal frequency (actual or estimate) [Hz]
    #           -verbose:    if 1, display information; else, don't. default: 1
    #           -plotflag:   If 1, plot; else, don't. default: 1
    #           -searchflag: If 0, iterative search skipped, just fit.
    #                        A quick choice you may want to try...
    #                        By default searchflag = 1.
    #
    # -The time vector is not necessarily linearly spaced (hence required here)
    # -If f is a guess, it has to be near true input signal frequency
    #
    # OUTPUTS:  -params:  A vector containing:
    #                     [offset amplitude freq(Hz) angle(rad)].
    #                     If not converged, params = [nan nan nan nan]
    #           -yest:    Estimated sinusoid:
    #                     offset+amplitude*cos(2*pi*frequency_Hz*t+angle_rad)
    #           -yres:    The residual
    #           -rmserr:  Root mean square of the residual
    #
    # For further information, consult IEEE Std 1057 and/or
    # IEEE Std 1241 documentation.
    #
    # type sinefit for demo
    # seed number for random generator. For repeatability...


    """
    if len(args) == 0   #demo
        clc
        fprintf('\ndemo1: creating sinusoid with\n\t- clock jitter,\n\t')
        fprintf('- phase noise,\n\t- additive noise and \n\t- harmonics\n ')
        N=pow2(12);          	%vector length
        fs = pi*1e3;            %sampling freq
        ts = 1/fs;              %sampling interval
        freq = (N/128-1)*fs/N;  %freq (Hz)
        phase = pi/2;          %phase (rad)
        offset = pi*2;          %offset (i.e mean)
        amplitude = pi/3;       %amplitude
        t = (0:(N-1))*ts;       %time vector
        std_jitter = 1e-2;  %standard deviation of jitter noise
        std_addnoi = 1e-1;  %standard deviation of  noise added to signal
        std_phase  = 1e-2;   %standard deviation of phase noise
        if oct_flag
            noise = randn(1,N);
        else
            noise = randn(s,1,N);
        end

        std1_noise = noise/std(noise);  % random vector with stdev = 1
        jit_noise = std_jitter*std1_noise;
        phase_noise = std_phase*std1_noise;
        add_noise = std_addnoi*std1_noise;
        w=2*pi*freq;
        t = t + ts*jit_noise;                          % add clock jitter
        A2 = amplitude*0.01;     % 2. harmonic ampl
        A3 = amplitude*0.02;     % 3. harmonic ampl
        yin = cos(w*t+phase+phase_noise);  % sinusoid with phase noise
        %add offset, noise & harmonics
        yin = offset+amplitude*yin+A2*yin.*yin+A3*yin.*yin.*yin+add_noise;
        figure
        params = sinefit(yin,t,freq,1,1);
        %params = sinefit(yin,t,freq,1,1,0);  %quick mode
        fprintf('\n\t\t\tpure sinusoid\tnoisy sinusoid fit')
        fprintf('\nOffset:   \t%3.4g\t\t\t%3.4g',offset,params(1))
        fprintf('\nAmpl:   \t%3.4g\t\t\t%3.4g',amplitude,params(2))
        fprintf('\nFreq:   \t%3.4g\t\t\t%3.4g [Hz]',freq,params(3))
        fprintf('\nPhase:\t\t%3.4g*pi\t\t%3.4g*pi [rad]\n',phase/pi,...
            params(4)/pi)
        fprintf('\nend demo1\n\n Press space for demo2')
        pause
        noiseq = randn(s,1,N);
        std1_noiseq = noise/std(noiseq);  % random vector with stdev = 1
        add_noise = 1e-2*(std1_noise+1i*std1_noiseq);
        fprintf('\n\ndemo2: phasor with noise and offset')
        yin = exp(1i*(w*t+phase));  % phasor with phase noise
        offset = offset + 1i*offset;
        yin = offset+amplitude*yin+add_noise;
        params = sinefit(yin,t,freq,1,1);
        fprintf('\n\t\t\tpure sinusoid\tnoisy sinusoid fit')
        fprintf('\nOffset:   \t%3.4g+j%3.4g\t%3.4g+j%3.4g',...
            real(offset),imag(offset),real(params(1)),imag(params(1)))
        fprintf('\nAmpl:   \t%3.4g\t\t\t%3.4g',amplitude,params(2))
        fprintf('\nFreq:   \t%3.4g\t\t\t%3.4g [Hz]',freq,params(3))
        fprintf('\nPhase:\t\t%3.4g*pi\t\t%3.4g*pi [rad]\n',phase/pi,...
            params(4)/pi)
        fprintf('\nend demo2\n')
        %clear params yest
        return
        end   %end demo
        """
    # convergence related parameters you can tweak:
    TOL = 1.e-15  # Normalized initial tolerance
    MTOL = 2.  # TOL relaxation multiplicand, MTOL > 1
    MAX_FUN = 22.  # MAX number of function calls
    MAX_ITER = 22.  # MAX number of iterations in one function call

    oct_flag = 0
    try:
        np.random.seed(1684)
    except:
        oct_flag = 1
        print ('Cannot create seed number')

    # varargin
    if not (oct_flag):
        if len(args) < 3 and len(args) > 0:
            print 'ERROR: at least three input params needed'
            return
        else:
            yin = args[0]
            t = args[1]
            f = args[2]

            verbose = 1
            plotflag = 0
            searchflag = 1

        if len(args) == 4:
            verbose = args[3]
            plotflag = 0
            searchflag = 1

        elif len(args) == 5:
            verbose = args[3]
            plotflag = args[4]
            searchflag = 1

        elif len(args) == 6:
            verbose = args[3]
            plotflag = args[4]
            searchflag = args[5]

        else:
            print('ERROR: too many input arguments')
            return

    # Convergence related stuff: Vector x0 will be created: x0=[A B C dw],
    # where A & B are sin & cos multiplicands, C is the offset and dw is
    # the angular frequency increment in iterative search.
    #
    # x0 is, of course, normalized in convergence calculation
    #(Fit for 1*yin converges as fast as for 100*yin)
    # ****  if max(abs(x0(i)-x0(i-1)))<TOL, the fit is complete.
    # ****  if not, multiply TOL with MTOL and retry maximum of MAX_FUN times.

    #plotting related
    if plotflag:
        N_per = 8  # number of periods plotted
        N_hist = 11  # number of bins in histogram

    if verbose:
        tim = time.time()

    N = len(yin)
    onevec = np.ones(N)

    # t does not have to be linearly spaced, so to estimate sampling time ts
    ts = rms(np.diff(t), N - 1)  # ts needed in freq normalization of converg. calc

    if MTOL < 0:
        print('ERROR: MTOL is a positive number > 1')

    elif MTOL < 1:
        MTOL = 1 / MTOL
        print('warning: MTOL should be > 1,')
        print('swiching to inverse value. XTOL = %f' % MTOL)

    # convergence related normalization
    rmsyin = rms(yin, N)
    TOL = rmsyin * TOL

    if len(args) < 6:
        searchflag = 1

    # here we go
    if not (searchflag):
        x0 = sinefit3par(yin, 2 * np.pi * f * t, onevec)
        x0 = np.append(x0, 2 * np.pi * f)  # not searching for freq
        success = 1
        if verbose:
            print('Quick mode: using 3-parameter sinefit')

    else:
        [x0, iter] = sinefit4par(yin, t, ts, 2 * np.pi * f, onevec, TOL, MAX_ITER)
        iter_total = iter
        iter_sinefit4par = 1  # first function call

        # success?
        if iter <= MAX_ITER:
            success = 1
            if verbose:
                print('\n\tConverged after ' + str(iter) + ' iterations\n')
        else:
            if verbose:
                print('\n\tincreasing TOL...')
            while iter > MAX_ITER and iter_sinefit4par <= MAX_FUN:
                if verbose:
                    print('.')
                # while number of function calls is < MAX_FUN, do:
                TOL = TOL * MTOL;  #increase TOL
                if oct_flag:
                    f = f + f / 80. * np.random.randn()
                else:
                    # reset(s)
                    f = f + f / 80. * np.random.randn()

                [x0, iter] = sinefit4par(yin, t, ts, 2 * np.pi * f, onevec, TOL, MAX_ITER)
                iter_total = iter_total + iter
                iter_sinefit4par = iter_sinefit4par + 1

            if iter > MAX_ITER:
                if verbose:
                    print('\nFailed to fit. The reasons could be:\n\t1. Bad '
                          'initial guess for frequency OR\n\t2. '
                          'the amplitude level is way below RMS noise floor OR'
                          '\n\t3. the fundamental frequency drifts (retry with a portion of input signal) OR'
                          '\n\t4. too small MAX_FUN, FUN_ITER, MTOL or TOL.\n\n')

                    print('%g function calls made, %g iterations allowed '
                          'in each.\n\n' % (iter_sinefit4par - 1, MAX_ITER))

                success = 0
                # return
            else:
                success = 1
                if verbose:
                    print('converged!\n')
                    print('\t%g function calls made,\n\t%g iterations allowed '
                          'in each.\n\n' % (iter_sinefit4par, MAX_ITER))

    # prep the output parameters
    A0 = x0[0]
    B0 = x0[1]
    C0 = x0[2]
    w = x0[3]
    sinedata = x0[0] * np.cos(x0[3] * t) + x0[1] * np.sin(x0[3] * t);
    yest = x0[2] + sinedata;

    if np.all(np.isreal(yin)):
        f_est = w / (2 * np.pi)
        A = np.sqrt(A0 * A0 + B0 * B0)
        phi = np.arctan(-B0 / A0)
        if A0 < 0:
            phi = phi + np.pi

        params = [C0, A, f_est, phi]
    else:
        f_est = np.real(w / (2 * np.pi))
        phi = np.angle(A0)
        # phi = angle(x0(1)/2+x0(2)/2)+pi/2;
        # temp = -B0/A0;
        # phi=atan(real(temp))+atan(imag(temp))-pi/4;
        # if real(A0)<0
        #     phi=phi+pi;
        # end
        params = [C0, abs(A0 / 2.) + abs(B0 / 2.), f_est, phi]

    yres = yin - yest
    rmserr = rms(yres, N)

    if verbose:
        t_elapsed = time.time() - tim;
        if t_elapsed < 60:
            print('\tTime elapsed: %g seconds\n' % t_elapsed)
        else:
            # this is not likely to happen
            print('\tTime elapsed: %g minutes and %g seconds\n' % np.floor(t_elapsed / 60), np.rem(t_elapsed, 60))

        if plotflag:
            print('\nplotting...')


    #plot or not
    #if plotflag == 1
    #    if isreal(yin)
    #        plotsinefit(N_hist,N_per,t,ts,yin,yest,yres,f_est,N,verbose)
    #    else
    #        figure
    #        plot(yin,'b')
    #        hold on
    #        plot(yest,'r')
    #        hold off
    #        xlabel('real')
    #        ylabel('imag')
    #        legend('data','fit',1)

    if not (success):
        params = ['nan', 'nan', 'nan', 'nan']
        yest = 'nan';
        yres = 'nan';
        rmserr = 'nan'


    # varargout
    # varargout[0]=params
    # varargout[1]=yest
    # varargout[2]=yres
    # varargout[3]=rmserr

    return [params, yest, yres, rmserr]


def sinefit3par(yin, wt, onevec):
    """3-parameter fit is used to create an initiall guess in sinefit4par
    """
    cosvec = np.cos(wt)
    sinvec = np.sin(wt)
    D0 = np.zeros((len(cosvec), 3))
    D0[:, 0] = cosvec
    D0[:, 1] = sinvec
    D0[:, 2] = onevec
    # x0=inv(D0.'*D0)*(D0.'*yin);
    if np.all(np.isreal(yin)):
        [Q, R] = np.linalg.qr(D0)
        x0 = np.linalg.solve(R, (np.dot(Q.T, yin)))  ### R\(Q.'*yin);
    else:
        x0 = np.linalg.lstsq(D0, yin)

    return x0


def sinefit4par(yin, t, ts, w, onevec, TOL, MAX_ITER):
    x0 = sinefit3par(yin, w * t, onevec)
    x0 = np.append(x0, 0)
    iter = 0
    success = 0
    while success == 0:
        w = w + x0[3]
        wt = w * t
        cosvec = np.cos(wt)
        sinvec = np.sin(wt)
        D0 = np.zeros((len(cosvec), 4))
        D0[:, 0] = cosvec
        D0[:, 1] = sinvec
        D0[:, 2] = onevec
        D0[:, 3] = -x0[0] * t * sinvec + x0[1] * t * cosvec
        x0old = x0

        # x0=inv(D0.'*D0)*(D0.'*yin);
        if np.all(np.isreal(yin)):
            [Q, R] = np.linalg.qr(D0)
            x0 = np.linalg.solve(R, np.dot(Q.T, yin))
        else:
            x0 = np.linalg.lstsq(D0, yin)

        iter += 1

        # error term with dw normalized
        temp = abs(x0 - x0old) * [1, 1, 1, ts]
        err = max(temp);

        if err < TOL or iter > MAX_ITER:  # if iter>MAX_ITER, increase TOL and
            success = 1  # re-visit this function later.

    x0[-1] = w  # place w in the position if w's increment

    return [x0, iter]


# function plotsinefit(N_hist,N_per,t,ts,yin,yest,yres,f,N,verbose)
#    [Nh,X] = hist([yin,yest],N_hist);
#    subplot(232)
#    barh(X,Nh)
#    title('Histograms')
#    legend('data','fit',1)
#    axis tight,xlabel('samples')
#
#    [Nh,X] = hist(yres,N_hist);
#    subplot(235)
#    barh(X,Nh,'k')
#    title('Residual histogram')
#    xlabel('samples')
#    axis tight
#
#    samples_per_period = abs(1/(f*ts));
#
#    new_N=ceil(N_per*samples_per_period);
#    if N >= new_N %if at least N_per (16) periods are found
#        N = new_N;
#        %selected_samples=1:N;
#        %t = t(selected_samples);
#        %yin = yin(selected_samples);
#        %yest = yest(selected_samples);
#        %yres = yres(selected_samples);
#    else
#        N_per=floor(N/samples_per_period);
#    end
#
#    selected_samples=1:N;
#    subplot(231)
#    plot(t(selected_samples),yin(selected_samples),'b',t(selected_samples),yest(selected_samples),'r')
#    axis tight
#    title([int2str(N_per),' periods plotted'])
#    legend('data','fit',1)
#    xlabel('t (s)')
#    ylabel('amplitude')
#    ylimy = get(gca,'YLim');
#
#    subplot(234)
#    plot(t(selected_samples),yres(selected_samples),'k')
#    xlabel('time (s)')
#    axis tight
#    title(['Residual, ',int2str(N_per),' periods'])
#    axis tight
#    xlabel('t (s)')
#    ylabel('amplitude')
#    ylime = get(gca,'YLim');
#
#    subplot(232)
#    set(gca,'YLim',ylimy)
#    subplot(235)
#    set(gca,'YLim',ylime)
#
#    tt=mod(t,1/f);
#    [tt,index]=sort(tt);
#    yy=yin(index);
#    yyest=yest(index);
#    yyresi = yres(index);
#    subplot(233)
#    plot(tt,yy,'b.',tt,yyest,'r','MarkerSize',4)
#    title('trace period plot')  % file ID: #22907
#    axis tight,legend('data','fit',1),xlabel('time, 1 period')
#    set(gca,'YLim',ylimy)
#    subplot(236)
#    plot(tt,yyresi,'k.','MarkerSize',4)
#    title('trace period plot, residual'), axis tight
#    xlabel('time, 1 period')
#    set(gca,'YLim',ylime)
#    if verbose
#        fprintf('done\n')
#    end

def rms(x, N):
    y = np.linalg.norm(x) / np.sqrt(N)
    return y


#t = np.linspace(0,10,1001)
#y = 5.0*np.cos(t*2*np.pi - np.pi/4)
#[params, yest, yres, rmserr] = sinefit(y,t,1.2,1,0,1)
#print params
#print yest
#print yres
#print rmserr
