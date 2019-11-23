import numpy as np
import cv2
import torch
import torchvision
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

# DO NOT CHANGE THE NAME OF THIS METHOD OR ITS INPUT OUTPUT BEHAVIOR

# INPUT CONVENTION
# filenames: a list of strings containing filenames of images

# OUTPUT CONVENTION
# The method must return a numpy array (not numpy matrix or scipy matrix) and a list of strings.
# Make sure that the length of the array and the list is the same as the number of filenames that 
# were given. The evaluation code may give unexpected results if this convention is not followed.

# https://nextjournal.com/gkoehler/pytorch-mnist
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = nn.Conv2d(1, 13, kernel_size=5)
        self.conv3 = nn.Conv2d(13, 26, kernel_size=3)
        self.conv2 = nn.Conv2d(26, 26, kernel_size=5)
        self.conv2_drop = nn.Dropout2d()
        self.fc1 = nn.Linear(650, 100)
        # self.fc2 = nn.Linear(150, 100)
        self.fc3 = nn.Linear(100, 26)

    def forward(self, x):
        x = F.relu(F.max_pool2d(self.conv1(x), 2))
        x = F.relu(F.max_pool2d(self.conv3(x), 4))
        x = F.relu(F.max_pool2d(self.conv2_drop(self.conv2(x)), 2))
        x = x.view(-1, 650)
        x = F.relu(self.fc1(x))
        x = F.dropout(x, training=self.training)
        x = self.fc3(x)
        return F.log_softmax(x, dim = 1)

def horizontal_sym(length, resize):
	size = (120 - length)//2
	temp = np.ones((150, size))
	resize = np.hstack((resize,temp*255))
	if length % 2 == 1:
		temp = np.hstack((temp, np.ones((150, 1))))
	resize = np.hstack((temp*255, resize))
	return resize

def vertical_sym(length, resize):
	size = (120 - length)//2
	temp = np.ones((size, 120))
	resize = np.vstack((resize,temp*255))
	if length % 2 == 1:
		temp = np.vstack((temp, np.ones((1, 120))))
	resize = np.vstack((temp*255, resize))
	return resize

def decaptcha( filenames ):
	l = len(filenames)
	numChars = np.array( [ 0 for i in range(l) ] )
	codes = []

	for index in range(l):
		image = cv2.imread(filenames[index])
		hsv_image = cv2.cvtColor(image,cv2.COLOR_BGR2HSV)
		_, gray_image = cv2.threshold(hsv_image[:,:,2], 203, 255, cv2.THRESH_BINARY)
		gray_image = np.asarray(gray_image)

		val = np.unique(np.where(gray_image==0)[1])
		length_val = len(val)
		start, end = val[0], -1
		segmented = []

		for i in range(length_val):
			if i == length_val-1:
				end = val[i]
			else:
				if val[i+1]-val[i] != 1:
					end = val[i]
			if end != -1:
				if end - start >= 20:
					resize = gray_image[:, start:end]
					resize = horizontal_sym(end-start, resize)

					hor_zeros = np.unique(np.where(resize==0)[0])
					resize = resize[hor_zeros[0]:hor_zeros[-1]+1, :]			
					resize = vertical_sym(hor_zeros[-1]-hor_zeros[0]+1, resize)
					resize = resize[np.newaxis, :, :]
					segmented.append(resize)

				if i < length_val-1:
					start = val[i+1]
				end = -1

		numChars[index] = len(segmented)

		# testing
		segmented = np.asarray(segmented)
		test_loader = torch.tensor(segmented)
		test_loader = test_loader.type(torch.FloatTensor)
		network = Net()
		network.load_state_dict(torch.load('model.pth'))
		network.eval()

		with torch.no_grad():
			data = test_loader
			output = network(data)
			pred = output.data.max(1, keepdim=True)[1]
			pred = pred.numpy()[:]
			pred = pred.tolist()
			pred = [chr(j+ord('A')) for sub in pred for j in sub]
			string = ""
			codes.append(string.join(pred))

	return (numChars, codes)