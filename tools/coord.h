//---------------------------------------------------------------------------

#ifndef coordH
#define coordH

#include <math.h>

namespace utilitis {

class TCoord {
    public:
        double x,y;
        TCoord(double X=0,double Y=0) {
            x=X,y=Y;
        }

        TCoord &operator=(const TCoord &c) {
            x =c.x;
            y=c.y;
            return *this;
        }

};

class CCoordException {

    public:
        int err;
        CCoordException(int em){ err = em; }
};

class TCoordChange {
    private:
        double xmin,xmax, gxmin, gxmax, ymin, ymax, gymin, gymax;
        double gxd, gyd, xd, yd;
    protected:
        void setXMin(double XMin) {
            if(XMin==xmax) {
                throw CCoordException(1);
            }
            xmin = XMin;
            if(xmax<xmin) {double t=xmax;xmax=xmin,xmin=t;}
            xd = 1.0/(xmax-xmin);
        }

        void setXMax(double XMax) {
            if(XMax==xmin) {
                throw CCoordException(1);
            }
            xmax = XMax;
            if(xmax<xmin) {double t=xmax;xmax=xmin,xmin=t;}
            xd = 1.0/(xmax-xmin);
        }

        void setgXMin(double XMin) {
            if(XMin==gxmax) {
                throw CCoordException(1);
            }
            gxmin = XMin;
            if(gxmax<gxmin) {double t=gxmax;gxmax=gxmin,gxmin=t;}
            gxd = 1.0/(gxmax-gxmin);
        }

        void setgXMax(double XMax) {
            if(XMax==gxmin) {
                throw CCoordException(1);
            }
            gxmax = XMax;
            if(gxmax<gxmin) {double t=gxmax;gxmax=gxmin,gxmin=t;}
            gxd = 1.0/(gxmax-gxmin);
        }

        void setYMin(double YMin) {
            if(YMin==ymax) {
                throw CCoordException(1);
            }
            ymin = YMin;
            if(ymax<ymin) {double t=ymax;ymax=ymin,ymin=t;}
            yd = 1.0/(ymax-ymin);
        }

        void setYMax(double YMax) {
            if(YMax==ymin) {
                throw CCoordException(1);
            }
            ymax = YMax;
            if(ymax<ymin) {double t=ymax;ymax=ymin,ymin=t;}
            yd = 1.0/(ymax-ymin);
        }

        void setgYMin(double YMin) {
            if(YMin==gymax) {
                throw CCoordException(1);
            }
            gymin = YMin;
            if(gymax<gymin) {double t=gymax;gymax=gymin,gymin=t;}
            gyd = 1.0/(gymax-gymin);
        }

        void setgYMax(double YMax) {
            if(YMax==gymin) {
                throw CCoordException(1);
            }
            gymax = YMax;
            if(gymax<gymin) {double t=gymax;gymax=gymin,gymin=t;}
            gyd = 1.0/(gymax-gymin);
        }

    public:
        TCoordChange(double Xmin=-1.0, double Xmax=1.0, double gXmin=0.0, double gXmax=1.0,
                     double Ymin=-1.0, double Ymax=1.0, double gYmin=0.0, double gYmax=1.0){

            if(Xmin==Xmax || gXmin==gXmax || Ymin==Ymax || gYmin==gYmax) {
                throw CCoordException(1);
            }

            double t;
            xmin = Xmin;
            xmax = Xmax;
            if(xmax<xmin) {t=xmax;xmax=xmin,xmin=t;}

            gxmin = gXmin;
            gxmax = gXmax;
            if(gxmax<gxmin) {t=gxmax;gxmax=gxmin,gxmin=t;}

            ymin = Ymin;
            ymax = Ymax;
            if(ymax<ymin) {t=ymax;ymax=ymin,ymin=t;}

            gymin = gYmin;
            gymax = gYmax;
            if(gymax<gymin) {t=gymax;gymax=gymin,gymin=t;}
            xd = 1.0/(xmax-xmin);
            gxd = 1.0/(gxmax-gxmin);
            yd = 1.0/(ymax-ymin);
            gyd = 1.0/(gymax-gymin);
        }

        TCoordChange(const TCoord &tl, const TCoord &br, const TCoord &gtl, const TCoord &gbr){
            if(tl.x==br.x || gtl.x==gbr.x || tl.y==br.y || gtl.y==gbr.y) {
                throw CCoordException(1);
            }

            double t;
            xmin = tl.x;
            xmax = br.x;
            if(xmax<xmin) {t=xmax;xmax=xmin,xmin=t;}

            gxmin = gtl.x;
            gxmax = gbr.x;
            if(gxmax<gxmin) {t=gxmax;gxmax=gxmin,gxmin=t;}

            ymin = tl.y;
            ymax = br.y;
            if(ymax<ymin) {t=ymax;ymax=ymin,ymin=t;}

            gymin = gtl.y;
            gymax = gbr.y;
            if(gymax<gymin) {t=gymax;gymax=gymin,gymin=t;}
            xd = 1.0/(xmax-xmin);
            gxd = 1.0/(gxmax-gxmin);
            yd = 1.0/(ymax-ymin);
            gyd = 1.0/(gymax-gymin);
        }

        double getX(double x) { // return in G coord
            return(gxmax-gxmin)*(x-xmin)*xd;
        }

        double getYi(double y) { // return in G coord
            return gymax-(gymax-gymin)*(y-ymin)*yd;
        }

        double getY(double y) { // return in G coord
            return(gymax-gymin)*(y-ymin)*yd;
        }

        double rgetX(double x) { // return in original coord
            return (x-gxmin)*(xmax -xmin)*gxd;
        }

        double rgetY(double y) { // return in original coord
            return (y-gymin)*(ymax - ymin)*gyd;
        }

        TCoord getXY(const TCoord &c) {
            return TCoord(getX(c.x), getY(c.y));
        }
        TCoord rgetXY(const TCoord &c) {
            return TCoord(rgetX(c.x), rgetY(c.y));
        }

        double getXMin() { return xmin; }

        double getXMax() { return xmax; }

        double getYMin() { return pow(10,ymin); }

        double getYMax() { return pow(10,ymax); }


        double getgXMin() { return gxmin; }

        double getgXMax() { return gxmax; }

        double getgYMin() { return gymin; }

        double getgYMax() { return gymax; }

};


class TCoordChangeLogY : public TCoordChange {
    public:
        TCoordChangeLogY(double Xmin=-1.0, double Xmax=1.0, double gXmin=0.0, double gXmax=1.0,
                     double Ymin=1e-6, double Ymax=10.0, double gYmin=0.0, double gYmax=1.0):TCoordChange(Xmin, Xmax, gXmin, gXmax, Ymin, Ymax, gYmin, gYmax){

            if(Ymin<=0 || Ymax<=0) {
                throw CCoordException(2);
            }

            setYMin(log10(Ymin));
            setYMax(log10(Ymax));
        }

        TCoordChangeLogY(const TCoord &tl, const TCoord &br, const TCoord &gtl, const TCoord &gbr):TCoordChange(tl,br,gtl,gbr){

            if(tl.y<=0 || br.y<=0) {
                throw CCoordException(2);
            }

            setYMin(log10(tl.y));
            setYMax(log10(br.y));
        }

        double getYi(double y) { // return in G coord
            return TCoordChange::getYi(log10(y));
        }

        double getY(double y) { // return in G coord
            return TCoordChange::getY(log10(y));
        }

        double rgetY(double y) { // return in original coord
            return TCoordChange::rgetY(log10(y));
        }

        TCoord getXY(const TCoord &c) {
            return TCoord(getX(c.x), getY(c.y));
        }
        TCoord rgetXY(const TCoord &c) {
            return TCoord(rgetX(c.x), rgetY(c.y));
        }

};

}
//---------------------------------------------------------------------------
#endif
