#ifndef GPMC_UI_H
#define GPMC_UI_H
/*
*
* file  gPMC_UI.cpp
* UI program for creating, inspecting and visualizing goPMC from commands.
*
* author Andreas Gravgaard Andersen
*
* last update on 26/7/2017
*
*/

// Qt stuff
#include <QtWidgets/QMainWindow>
#include "ui_gPMC.h"
// ITK stuff
#include "itkImage.h"
// My stuff
#include "YK16GrayImage.h"
#include "AG17RGBAImage.h"

#include "scatterdatamodifier.h"
#include <QtDataVisualization/q3dscatter.h>
#include <QtGui/QScreen>

#define DEFAULT_LABEL_SIZE1 512
#define DEFAULT_LABEL_SIZE2 256
#define DEFAULT_LABEL_SIZE3 256

class qyklabel;

enum enViewArrange{
	AXIAL_FRONTAL_SAGITTAL = 0,
	FRONTAL_SAGITTAL_AXIAL,
	SAGITTAL_AXIAL_FRONTAL,
};

enum enPLANE{
	PLANE_AXIAL = 0,
	PLANE_FRONTAL,
	PLANE_SAGITTAL,
};

typedef itk::Image<unsigned short, 3> UShortImageType;
typedef itk::Image<short, 3> ShortImageType;
typedef itk::Image<float, 3> FloatImageType;

class gPMC_UI : public QMainWindow
{
	Q_OBJECT
		;

public:
	gPMC_UI(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	~gPMC_UI();

	void whenFixedImgLoaded();
	void MousePressedRight(int wndIdx, qyklabel* pWnd);
	void UpdateSplit(int viewIdx, qyklabel* pOverlapWnd);
	void initOverlapWndSize();
	void shiftSliceSlider();
	void updateSliceLabel();
	void Draw2DFrom3DDouble(UShortImageType::Pointer& spFixedImg, UShortImageType::Pointer& spMovingImg, enPLANE direction, double pos, YK16GrayImage& YKFixed, YK16GrayImage& YKMoving);
	void Draw2DFrom3DDouble(UShortImageType::Pointer& spFixedImg, UShortImageType::Pointer& spMovingImg, enPLANE direction, double pos, AG17RGBAImage& YKFixed, AG17RGBAImage& YKMoving);

	void childEvent(QChildEvent *event);
	bool eventFilter(QObject *target, QEvent *event);
	void UpdateListOfComboBox(int idx);
	void LoadImgFromComboBox(int idx, QString& strSelectedComboTxt);
	void load_dcm_into_ScatterWidget();

	public slots:
	void SLT_CrntPosGo();
	void SLT_DrawImageWhenSliceChange(); //upper level drawing: big calculation
	void SLT_DrawImageInFixedSlice();//lower level Drawing func.

	void SLT_UpdateSplit1();//lower level Drawing func. //Mouse Move even
	void SLT_UpdateSplit2();//lower level Drawing func.//Mouse Move even
	void SLT_UpdateSplit3();//lower level Drawing func.//Mouse Move even

	void SLT_CancelMouseAction();

	void SLT_MouseWheelUpdate1();//wheel
	void SLT_MouseWheelUpdate2();
	void SLT_MouseWheelUpdate3();
	void SLT_MousePressedLeft1();//Press right
	void SLT_MousePressedLeft2();
	void SLT_MousePressedLeft3();
	void SLT_MousePressedRight1();//Press right
	void SLT_MousePressedRight2();
	void SLT_MousePressedRight3();
	void SLT_MouseReleasedLeft1();//Rel left
	void SLT_MouseReleasedLeft2();
	void SLT_MouseReleasedLeft3();
	void SLT_MouseReleasedRight1();//Rel rigth
	void SLT_MouseReleasedRight2();
	void SLT_MouseReleasedRight3();
	void SLT_ChangeView();//3 toggle button
	void SLT_KeyMoving(bool bChecked);
	void SLT_BringFocusToEnableArrow(bool bChecked);

	void SLT_recalcDose();
	void SLT_pltBeamGeo();
	//Image selection event from Combobox
	void SLT_FixedImageSelected(QString selText); //here, when fixed_image_loaded function will be called
	void SLT_MovingImageSelected(QString selText);//here, when mvoing_image_loaded function will be called
	void SLT_RestoreImageSingle();
	void SLT_RestoreImageAll();
	void SLTM_LoadDICOMdir();

private:

public:
	YK16GrayImage m_YKImgFixed[3]; //CBCT in this study
	YK16GrayImage m_YKImgMoving[3]; //CBCT in this study
	YK16GrayImage m_YKDisp[3]; //CBCT in this study
	AG17RGBAImage m_DoseImgFixed[3]; //CBCT in this study
	AG17RGBAImage m_DoseImgMoving[3]; //CBCT in this study
	AG17RGBAImage m_AGDisp_Overlay[3]; //CBCT in this study
	bool dose_loaded = false;
	ScatterDataModifier *m_modifier;

	UShortImageType::Pointer m_spFixed;//pointer only, for display
	UShortImageType::Pointer m_spMoving;//pointer only, for display

	UShortImageType::Pointer m_spFixedDose;//pointer only, for display
	UShortImageType::Pointer m_spMovingDose;//pointer only, for display

	UShortImageType::Pointer m_spRawReconImg;
	UShortImageType::Pointer m_spRefCTImg;

	int m_enViewArrange;
	QString m_strPathDirDefault = QString("./");

	bool m_bPressedLeft[3];//Left Mouse Pressed but not released
	bool m_bPressedRight[3];
	QPoint m_ptWindowLevelStart;//data point

	QPoint m_ptPanStart;//data point

	QPoint m_ptTmpOriginalDataOffset;
	int m_iTmpOriginalW;
	int m_iTmpOriginalL;

	Ui::gPMC ui;
};

#endif