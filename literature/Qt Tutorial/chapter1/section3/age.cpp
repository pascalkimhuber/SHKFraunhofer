#include <QApplication>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QSlider>
#include <QSpinBox>
#include <QPushButton>

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	
	QWidget *window = new QWidget;
	QWidget *func = new QWidget(window);
	window->setWindowTitle("Enter Your Age");

	QSpinBox *spinBox = new QSpinBox;
	QSlider *slider = new QSlider(Qt::Horizontal);
	QPushButton *button = new QPushButton("Cancel");
	
	spinBox->setRange(0, 130);
	slider->setRange(0, 130);

	QObject::connect(spinBox, SIGNAL(valueChanged(int)), slider, SLOT(setValue(int)));
	QObject::connect(slider, SIGNAL(valueChanged(int)), spinBox, SLOT(setValue(int)));
	QObject::connect(button, SIGNAL(clicked()), &app, SLOT(quit()));
	spinBox->setValue(35);

	QHBoxLayout *layout = new QHBoxLayout;
	QVBoxLayout *globalLayout = new QVBoxLayout;
	layout->addWidget(spinBox);
	layout->addWidget(slider);
	func->setLayout(layout);
	globalLayout->addWidget(func);
	globalLayout->addWidget(button);
	window->setLayout(globalLayout);

	window->show();

	return app.exec();
}
